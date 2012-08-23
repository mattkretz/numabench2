/*{{{
    Copyright (C) 2010-2012 Matthias Kretz <kretz@kde.org>

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
    02110-1301, USA.

}}}*/
#include <Vc/Vc>/*{{{*/
#include "benchmark.h"
#include <sys/mman.h>
#include <fstream>
#include <thread>
#include <mutex>
#include <atomic>
#include <Vc/cpuid.h>
#include "cpuset.h"
#ifndef NO_LIBNUMA
#include "numaif.h"
#include "numa.h"
#endif

using Vc::CpuId;
extern "C" void numa_init();
/*}}}*/
typedef float Scalar;
typedef Vc::Vector<Scalar> Vector;
typedef Scalar *__restrict__ Memory;
enum Constants {/*{{{*/
    PageSize = 4096,
    CacheLineSize = 64,

    ScalarsInPage = PageSize / sizeof(Scalar),
    VectorsInPage = PageSize / sizeof(Vector),
    ScalarsInCacheLine = CacheLineSize / sizeof(Scalar),
    VectorsInCacheLine = CacheLineSize / sizeof(Vector)
};
constexpr size_t GiB = 1024ull * 1024ull * 1024ull;
/*}}}*/
struct TestArguments/*{{{*/
{
    Memory mem;
    Timer *__restrict__ timer;
    size_t offset;
    size_t size;
    int repetitions;
};/*}}}*/
typedef void (*TestFunction)(const TestArguments &args);
struct MemoryRange/*{{{*/
{
    MemoryRange(Memory s, size_t size) : start(s), end(start + size) {}
    const Memory start, end;
};/*}}}*/
struct CpuRange/*{{{*/
{
    int first;
    int last;
    int step;
};/*}}}*/
inline std::ostream &operator<<(std::ostream &out, const CpuRange &range)/*{{{*/
{
    return out << "Range " << range.first << " - " << range.last << ", step " << range.step << '\n';
}
/*}}}*/
static int maxThreadCount()/*{{{*/
{
    cpu_set_t cpumask;
    sched_getaffinity(0, sizeof(cpu_set_t), &cpumask);
    return cpuCount(&cpumask);
}/*}}}*/

class OneWaitsForN/*{{{*/
{
private:
    std::mutex mutex;
    std::condition_variable wait;
    std::atomic<int> busyCount;
public:
    OneWaitsForN(int count)
        : busyCount(count)
    {
    }

    void oneReady()
    {
        if (busyCount-- == 1) {
            std::unique_lock<std::mutex> lock(mutex);
            wait.notify_one();
        }
    }

    void setBusy(int count)
    {
        busyCount = count;
    }

    void waitForAll()
    {
        if (busyCount > 0) {
            std::unique_lock<std::mutex> lock(mutex);
            if (busyCount > 0) {
                wait.wait(lock);
            }
        }
    }
};/*}}}*/
class ThreadData/*{{{*/
{
    std::mutex m_mutex;
    std::condition_variable_any &m_wait;
    OneWaitsForN &m_waitForEnd;
    std::atomic<bool> m_exit;
    bool m_disabled;
    int m_cpuId;
    Timer m_timer;
    std::thread m_thread;
    TestFunction m_testFunction;
    TestArguments m_arguments;

    public:
        ThreadData(OneWaitsForN *waitForEnd, std::condition_variable_any *wait) // called from main thread
            : m_wait(*wait),
            m_waitForEnd(*waitForEnd),
            m_exit(false),
            m_disabled(false),
            m_cpuId(-1),
            m_thread(ThreadData::callMainLoop, this)
        {}

        const Timer &timer() const { return m_timer; }

        void setPinning(int cpuid) // called from main thread
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_disabled = false;
            m_cpuId = cpuid;
        }

        void disable()
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_disabled = true;
        }
        void enable()
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_disabled = false;
        }
        bool isEnabled() const
        {
            return !m_disabled;
        }

        void exit() // called from main thread
        {
            m_exit = true;
        }

        void join() // called from main thread
        {
            m_thread.join();
        }

        void setTestFunction(TestFunction f)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_testFunction = f;
        }

        void setParameters(TestArguments args)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_arguments = args;
            m_arguments.timer = &m_timer;
        }

        static void callMainLoop(ThreadData *data) // thread
        {
            data->mainLoop();
        }

        void barrier()
        {
            // this function will block until/unless the mainLoop is in m_wait.wait
            std::lock_guard<std::mutex> lock(m_mutex);
        }

    private:
        void mainLoop() // thread
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            do {
                m_waitForEnd.oneReady();

                // wait for the signal to start
                m_wait.wait(m_mutex);

                if (m_exit) {
                    break;
                } else if (m_disabled) {
                    continue;
                } else if (m_cpuId >= 0) {
                    // first pin the thread to a single core/cpu
                    pinToCpu(m_cpuId);
                    m_cpuId = -1;
                    continue;
                }


                // do the work
                m_testFunction(m_arguments);
            } while (!m_exit);
            m_waitForEnd.oneReady();
        }
};/*}}}*/
class ThreadPool/*{{{*/
{
    OneWaitsForN m_waitForEnd;
    std::condition_variable_any m_waitForStart;
    std::vector<std::shared_ptr<ThreadData>> m_workers;

public:
    ThreadPool(int _size = maxThreadCount())
        : m_waitForEnd(_size),
        m_workers(_size)
    {
        for (int i = 0; i < _size; ++i) {
            m_workers[i] = std::make_shared<ThreadData>(&m_waitForEnd, &m_waitForStart);
        }
    }

    void waitReady()
    {
        m_waitForEnd.waitForAll();
    }

    void wakeAllWorkers()
    {
        m_waitForEnd.waitForAll();
        m_waitForEnd.setBusy(m_workers.size());
        for (auto &t : m_workers) {
            t->barrier();
        }
        m_waitForStart.notify_all();
    }

    /*void setPinning(int firstCpu)
    {
        for (size_t i = 0; i < m_workers.size(); ++i) {
            m_workers[i]->setPinning(firstCpu + i);
        }
    }
    void setPinning(std::vector<int> cpus)
    {
        assert(m_workers.size() == cpus.size());
        for (size_t i = 0; i < m_workers.size(); ++i) {
            m_workers[i]->setPinning(cpus[i]);
        }
    }*/
    void setPinning(CpuRange range)
    {
        unsigned int i = 0;
        for (int id = range.first; id <= range.last; id += range.step) {
            m_workers[i++]->setPinning(id);
        }
        for (; i < m_workers.size(); ++i) {
            m_workers[i]->disable();
        }
        wakeAllWorkers(); // have the threads call sched_setaffinity
    }

    void setTestFunction(TestFunction f)
    {
        for (auto &t : m_workers) {
            t->setTestFunction(f);
        }
    }

    template<typename OffsetFunction>
    void executeWith(Memory mem, OffsetFunction offset, size_t _size, int repetitions)
    {
        for (auto &t : m_workers) {
            t->setParameters({ mem, nullptr, offset(), _size, repetitions });
        }
        wakeAllWorkers();
    }

    template<typename F> void eachTimer(F f) const
    {
        for (const auto &t : m_workers) {
            if (t->isEnabled()) {
                f(t->timer());
            }
        }
    }

    ~ThreadPool()
    {
        for (auto &t : m_workers) {
            t->exit();
        }
        m_waitForStart.notify_all();
        for (auto &t : m_workers) {
            t->join();
        }
    }
};/*}}}*/

template<size_t SliceSize = 1> struct TestDefaults/*{{{*/
{
    /// in #Scalars
    static constexpr size_t sliceSizeT() { return SliceSize * GiB / sizeof(Scalar); }
    /// in Bytes
    static std::vector<size_t> sizes() { return { SliceSize * GiB, CpuId::L3Data() / 2, CpuId::L2Data() / 2, CpuId::L1Data() / 2 }; }
    /// in #Scalars
    static std::vector<size_t> offsetsPerThread() {
        return { 0,
            CpuId::L1Data() / sizeof(Scalar),
            sliceSizeT() / 128,
            sliceSizeT() / 128 + 7 * ScalarsInCacheLine
        };
    }
    static constexpr double interpretFactor() { return 1.; }
    static constexpr const char *interpretUnit() { return "Byte"; }
    /// in #Scalars
    static constexpr size_t stride() { return 1 * GiB / sizeof(Scalar); }
    static constexpr const char *description() { return "undocumented"; }
};/*}}}*/
template<bool Streaming> struct TestWriteBase : public TestDefaults<1>/*{{{*/
{
    static constexpr const char *name() { return "write"; }
    static void run(const TestArguments &args)
    {
        const size_t offset = args.offset + args.size > sliceSizeT() ? args.offset : 0;
        const MemoryRange mRange[2] = {
            { args.mem + args.offset, args.size - offset },
            { args.mem, offset }
        };
        args.timer->start();
        const Vector one = Vector::One();
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (int i = 0; i < 2; ++i) {
                for (Memory m = mRange[i].start + 3 * Vector::Size; m < mRange[i].end; m += 4 * Vector::Size) {
                    if (Streaming) {
                        one.store(m - 3 * Vector::Size, Vc::Streaming);
                        one.store(m - 2 * Vector::Size, Vc::Streaming);
                        one.store(m - 1 * Vector::Size, Vc::Streaming);
                        one.store(m - 0 * Vector::Size, Vc::Streaming);
                    } else {
                        one.store(m - 3 * Vector::Size);
                        one.store(m - 2 * Vector::Size);
                        one.store(m - 1 * Vector::Size);
                        one.store(m - 0 * Vector::Size);
                    }
                }
            }
        }
        args.timer->stop();
    }
};/*}}}*/
struct TestWrite : public TestWriteBase<false>/*{{{*/
{
    static constexpr const char *name() { return "write"; }
    static constexpr const char *description() { return
        "The write benchmark naïvely writes constant data over the given memory.\n\n"
        "Depending on the microarchitecture non-streaming writes require the memory\n"
        "controller to issue loads before the stores, in order to initialize the cache\n"
        "lines that will be overwritten.";
    }
};/*}}}*/
struct TestStream : public TestWriteBase<true>/*{{{*/
{
    static constexpr const char *name() { return "stream"; }
    static constexpr const char *description() { return
        "The stream benchmark writes constant data over the given memory using\n"
        "non-temporal stores.\n";
    }
};/*}}}*/
struct TestBzero : public TestDefaults<1>/*{{{*/
{
    static constexpr const char *name() { return "bzero"; }
    static void run(const TestArguments &args)
    {
        if (args.offset + args.size > sliceSizeT()) {
            args.timer->start();
            for (int rep = 0; rep < args.repetitions; ++rep) {
                // TODO exchange the two? because bzero initializes backwards AFAIK
                bzero(args.mem + args.offset, (args.size - args.offset) * sizeof(Scalar));
                bzero(args.mem, args.offset * sizeof(Scalar));
            }
            args.timer->stop();
        } else {
            args.timer->start();
            for (int rep = 0; rep < args.repetitions; ++rep) {
                bzero(args.mem + args.offset, args.size * sizeof(Scalar));
            }
            args.timer->stop();
        }
    }
    static constexpr const char *description() { return
        "The bzero benchmark simply calls bzero(3) for the given memory ranges.\n\n"
        "A typical bzero implementation uses memset, which in turn should be implemented\n"
        "to initialize memory with streaming stores, except for the first few kilobytes.\n"
        "Interpretation of the bzero results is always dependent on the specific\n"
        "implementation and thus should with great care.";
    }
};/*}}}*/
struct TestAddOneStrided : public TestDefaults<8>/*{{{*/
{
    // In the "big test" iterate over 8 GiB with strides of the size of the L2 cache. For an L2
    // cache of 2 MiB that results in just 4096 loads.
    // In the "L3 test" the loads should be serviced by the L3 cache.
    static std::vector<size_t> sizes() { return { 8 * GiB, CpuId::L3Data() }; }
    static constexpr const char *name() { return "add 1 w/ large strides"; }
    static constexpr double interpretFactor() { return 1./1024./sizeof(Scalar); }
    static constexpr const char *interpretUnit() { return "Add"; }
    static void run(const TestArguments &args)
    {
        // there are (args.size / Stride) many adds in the inner loop
        // the inner loop is repeated Stride/1024 times
        // => there will roughly be args.size adds in total
        Memory mStart = args.mem;
        Memory const mEnd = mStart + args.size;
        const size_t Stride = CpuId::L2Data() / sizeof(Scalar);
        if (Stride >= args.size) {
            const Vector one(1);
            const int repetitions = args.repetitions / 1024;
            args.timer->start();
            for (int rep = 0; rep < repetitions; ++rep) {
                for (Memory m = mStart + 3 * Vector::Size; m < mEnd; m += 4 * Vector::Size) {
                    const Vector tmp0 = Vector(m - 3 * Vector::Size) + one;
                    const Vector tmp1 = Vector(m - 2 * Vector::Size) + one;
                    const Vector tmp2 = Vector(m - 1 * Vector::Size) + one;
                    const Vector tmp3 = Vector(m - 0 * Vector::Size) + one;
                    tmp0.store(m - 3 * Vector::Size);
                    tmp1.store(m - 2 * Vector::Size);
                    tmp2.store(m - 1 * Vector::Size);
                    tmp3.store(m - 0 * Vector::Size);
                }
            }
            args.timer->stop();
        } else {
            const Scalar one(1);
            const int repetitions = args.repetitions * Stride / 1024;
            args.timer->start();
            for (int rep = 0; rep < repetitions; ++rep) {
                for (Memory m = mStart; m < mEnd; m += Stride) {
                    *m += one;
                }
                mStart += ScalarsInCacheLine;
            }
            args.timer->stop();
        }
    }
    static constexpr const char *description() { return
        "The add 1 w/ large strides benchmark tries to render the caches useless and thus\n"
        "limit the bandwidth to some kind of worst-case scenario.\n\n"
        "The benchmark modifies values in memory at strides of exactly the size of the L2\n"
        "cache. Thus the L1/L2 caches can only store as many values as the associativity\n"
        "of the cache allows. Typically this reduces the effective cache to the size of\n"
        "8/12/16 cache lines (512/768/1024 Bytes). Subsequent memory reads have to wait\n"
        "until the data in the cache has been evicted. Thus the latency of memory loads\n"
        "(and NUMA interconnects) will show greatly in this benchmark.";
    }
};/*}}}*/
struct TestAddOneStrided2 : public TestDefaults<8>/*{{{*/
{
    static std::vector<size_t> sizes() { return { 8 * GiB }; }
    static constexpr const char *name() { return "add 1 w/ optimized strides"; }
    static constexpr double interpretFactor() { return 1./sizeof(Scalar); }
    static constexpr const char *interpretUnit() { return "Add"; }
    static void run(const TestArguments &args)
    {
        const Vector one(1);
        Memory mStart = args.mem;
        Memory const mEnd = mStart + args.size;
        constexpr size_t Stride = (32 * 1024 + CacheLineSize) / sizeof(Scalar);
        Memory mPartEnd = mStart + Stride;
        const int repetitions = args.repetitions;
        args.timer->start();
        for (int rep = 0; rep < repetitions; ++rep) {
            while (mPartEnd <= mEnd) {
                for (Memory m = mStart + Vector::Size; m < mPartEnd; m += 2 * Vector::Size) {
                    const Vector tmp0 = Vector(m - 1 * Vector::Size) + one;
                    const Vector tmp1 = Vector(m - 0 * Vector::Size) + one;
                    const Vector tmp2 = Vector(m + Stride - 1 * Vector::Size) + one;
                    const Vector tmp3 = Vector(m + Stride - 0 * Vector::Size) + one;
                    tmp0.store(m - 1 * Vector::Size);
                    tmp1.store(m - 0 * Vector::Size);
                    tmp2.store(m + Stride - 1 * Vector::Size);
                    tmp3.store(m + Stride - 0 * Vector::Size);
                }
                mStart += 2 * Stride;
                mPartEnd += 2 * Stride;
            }
        }
        args.timer->stop();
    }
};/*}}}*/
template<bool Prefetch> struct TestAddOneBase : public TestDefaults<1>/*{{{*/
{
    static constexpr double interpretFactor() { return 2.; } // every memory access is both load and store, thus twice the bandwidth
    static void run(const TestArguments &args)
    {
        const Vector one = Vector::One();

        const size_t offset = args.offset + args.size > sliceSizeT() ? args.offset : 0;
        const MemoryRange mRange[2] = {
            { args.mem + args.offset, args.size - offset },
            { args.mem, offset }
        };

        args.timer->start();
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (int i = 0; i < 2; ++i) {
                for (Memory m = mRange[i].start + 3 * Vector::Size; m < mRange[i].end; m += 4 * Vector::Size) {
                    if (Prefetch) {
                        Vc::prefetchForModify(m + 4032 / sizeof(Scalar));
                    }
                    const Vector tmp0 = Vector(m - 3 * Vector::Size) + one;
                    const Vector tmp1 = Vector(m - 2 * Vector::Size) + one;
                    const Vector tmp2 = Vector(m - 1 * Vector::Size) + one;
                    const Vector tmp3 = Vector(m - 0 * Vector::Size) + one;
                    tmp0.store(m - 3 * Vector::Size);
                    tmp1.store(m - 2 * Vector::Size);
                    tmp2.store(m - 1 * Vector::Size);
                    tmp3.store(m - 0 * Vector::Size);
                }
            }
        }
        args.timer->stop();
    }
};/*}}}*/
struct TestAddOne : public TestAddOneBase<false>/*{{{*/
{
    static constexpr const char *name() { return "add 1"; }
};/*}}}*/
struct TestAddOnePrefetch : public TestAddOneBase<true>/*{{{*/
{
    static constexpr const char *name() { return "add 1 w/ prefetch"; }
};/*}}}*/
template<bool Prefetch> struct TestReadBase : public TestDefaults<1>/*{{{*/
{
    static void run(const TestArguments &args)
    {
        const size_t offset = args.offset + args.size > sliceSizeT() ? args.offset : 0;
        const MemoryRange mRange[2] = {
            { args.mem + args.offset, args.size - offset },
            { args.mem, offset }
        };
        args.timer->start();
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (int i = 0; i < 2; ++i) {
                for (Memory m = mRange[i].start + 7 * Vector::Size; m < mRange[i].end; m += 8 * Vector::Size) {
                    if (Prefetch) {
                        Vc::prefetchForOneRead(m + 4032 / sizeof(Scalar));
                    }
                    const Vector v0(m - 7 * Vector::Size);
                    const Vector v1(m - 6 * Vector::Size);
                    const Vector v2(m - 5 * Vector::Size);
                    const Vector v3(m - 4 * Vector::Size);
                    const Vector v4(m - 3 * Vector::Size);
                    const Vector v5(m - 2 * Vector::Size);
                    const Vector v6(m - 1 * Vector::Size);
                    const Vector v7(m - 0 * Vector::Size);
                    Vc::forceToRegisters(v0, v1, v2, v3, v4, v5, v6, v7);
                }
            }
        }
        args.timer->stop();
    }
};/*}}}*/
struct TestRead : public TestReadBase<false>/*{{{*/
{
    static constexpr const char *name() { return "read"; }
};/*}}}*/
struct TestReadPrefetch : public TestReadBase<true>/*{{{*/
{
    static constexpr const char *name() { return "read w/ prefetch"; }
};/*}}}*/
/** testReadLatency {{{
 * We want to measure the latency of a read from memory. To achieve this we read with a stride of
 * PageSize bytes. Then the hardware prefetcher will not do any prefetches and every load will hit a
 * cold cache line. To increase the working size the test then starts over but with an offset of one
 * cache line:
 * [x                               x                               ...]
 * [        x                               x                       ...]
 * [                x                               x               ...]
 * [                        x                               x       ...]
 */
struct TestReadLatency : public TestDefaults<1>
{
    static constexpr const char *name() { return "read latency"; }
    static constexpr double interpretFactor() { return 1. / VectorsInCacheLine; }
    static constexpr const char *interpretUnit() { return "Read"; }
    static void run(const TestArguments &args)
    {
        Memory const mStart = args.mem;
        Memory const mEnd = mStart + args.size;
        Memory const mPageEnd = mStart + ScalarsInPage;
        args.timer->start();
        if (((mEnd - mStart) / ScalarsInCacheLine) & 1) {
            for (int rep = 0; rep < args.repetitions; ++rep) {
                for (Memory mCacheLine = mStart; mCacheLine < mPageEnd; mCacheLine += ScalarsInCacheLine) {
                    for (Memory m = mCacheLine; m < mEnd; m += ScalarsInPage) {
                        //asm volatile("lfence");
                        asm volatile("" :: "d"(*m));
                    }
                }
            }
        } else if (((mEnd - mStart) / ScalarsInCacheLine) & 2) {
            for (int rep = 0; rep < args.repetitions; ++rep) {
                for (Memory mCacheLine = mStart; mCacheLine < mPageEnd; mCacheLine += ScalarsInCacheLine) {
                    for (Memory m = mCacheLine; m < mEnd - ScalarsInPage; m += 2 * ScalarsInPage) {
                        //asm volatile("lfence");
                        asm volatile("" :: "d"(*m));
                        asm volatile("" :: "d"(*(m + ScalarsInPage)));
                    }
                }
            }
        } else {
            for (int rep = 0; rep < args.repetitions; ++rep) {
                for (Memory mCacheLine = mStart; mCacheLine < mPageEnd; mCacheLine += ScalarsInCacheLine) {
                    for (Memory m = mCacheLine; m < mEnd - 3 * ScalarsInPage; m += 4 * ScalarsInPage) {
                        //asm volatile("lfence");
                        asm volatile("" :: "d"(*m));
                        asm volatile("" :: "d"(*(m + ScalarsInPage)));
                        asm volatile("" :: "d"(*(m + 2 * ScalarsInPage)));
                        asm volatile("" :: "d"(*(m + 3 * ScalarsInPage)));
                    }
                }
            }
        }
        args.timer->stop();
    }
};/*}}}*/

template<typename T> struct convertStringTo/*{{{*/
{
    explicit convertStringTo(const std::string &s);
    operator T() { return m_data; }
    T m_data;
};
/*}}}*/
template<> convertStringTo<int>::convertStringTo(const std::string &s) : m_data(atoi(s.c_str())) {}/*{{{*/
template<> convertStringTo<unsigned int>::convertStringTo(const std::string &s) : m_data(atoi(s.c_str())) {}
template<> convertStringTo<long>::convertStringTo(const std::string &s) : m_data(atol(s.c_str())) {}
template<> convertStringTo<unsigned long>::convertStringTo(const std::string &s) : m_data(atol(s.c_str())) {}
template<> convertStringTo<long long>::convertStringTo(const std::string &s) : m_data(atoll(s.c_str())) {}
template<> convertStringTo<unsigned long long>::convertStringTo(const std::string &s) : m_data(atoll(s.c_str())) {}
template<> convertStringTo<std::string>::convertStringTo(const std::string &s) : m_data(s) {}/*}}}*/
template<typename T> static T valueForArgument(const char *name, T defaultValue)/*{{{*/
{
    ArgumentVector::iterator it = std::find(g_arguments.begin(), g_arguments.end(), name);
    if (it != g_arguments.end()) {
        ++it;
        if (it != g_arguments.end()) {
            return convertStringTo<T>(*it);
        }
    }
    return defaultValue;
}/*}}}*/
std::vector<CpuRange> parseOnlyCpus()/*{{{*/
{
    enum State {
        ReadFirst, ReadLast, ReadStep
    };
    const std::string cpusStrings = valueForArgument("--cores", std::string());
    std::vector<CpuRange> r;
    CpuRange range = { 0, 0, 1 };
    State state = ReadFirst;
    const int rangeMax = maxThreadCount();
    for (const auto &c : cpusStrings) {
        if (c >= '0' && c <= '9') {
            switch (state) {
            case ReadFirst:
                range.first = range.first * 10 + (c - '0');
                break;
            case ReadLast:
                range.last = range.last * 10 + (c - '0');
                break;
            case ReadStep:
                range.step = range.step * 10 + (c - '0');
                break;
            }
        } else if (c == '-') {
            state = ReadLast;
        } else if (c == ':') {
            state = ReadStep;
            range.step = 0;
        } else if (c == ',') {
            if (state == ReadFirst) {
                range.last = range.first;
            } else {
                state = ReadFirst;
            }
            if (range.first < 0 || range.first >= rangeMax || range.last < range.first || range.last >= rangeMax) {
                std::cerr << "invalid cores: " << range.first << '-' << range.last << std::endl;
                std::exit(1);
            }
            r.push_back(range);
            range = { 0, 0, 1 };
        }
    }
    r.push_back(range);
    return r;
}
/*}}}*/
static size_t largestMemorySize()/*{{{*/
{
    using namespace std;
    fstream meminfo("/proc/meminfo", fstream::in);
    string tmp;
    size_t totalMem, freeMem;
    meminfo >> tmp >> totalMem >> tmp >> tmp >> freeMem;
    meminfo.close();
    return freeMem * 1024;
}/*}}}*/
SET_HELP_TEXT(/*{{{*/
        "  --listTests\n"
        "  --firstCpu <id>\n"
        "  --cpuStep <id>\n"
        "  --size <GiB>\n"
        "  --only <test function>\n"
        "  --cores <firstId-lastId[:step][,firstId-lastId[:step][...]]>\n"
        );/*}}}*/
static std::string prettyBytes(size_t bytes)/*{{{*/
{
    std::stringstream ss0;
    constexpr char const *Prefix[] = { "", "Ki", "Mi", "Gi", "Ti" };
    int prefixI = 0;
    if (bytes != 0) {
        while (bytes % 1024 == 0) {
            bytes /= 1024;
            ++prefixI;
        }
    }
    ss0 << bytes << ' ' << Prefix[prefixI] << "B";
    return ss0.str();
}/*}}}*/

class BenchmarkRunner/*{{{*/
{
private:
    const std::vector<CpuRange> m_coreIds;
    int m_threadCount;
    ThreadPool m_threadPool;
    const size_t m_maxMemorySize;
    const size_t m_memorySize;
    const std::string m_only;
    Memory m_memory;

    template<typename Test> void executeTest();
    void executeAllTests();

public:
    BenchmarkRunner();
};/*}}}*/
template<typename T> std::string toString(const T &x)/*{{{*/
{
    std::stringstream ss;
    ss << x;
    return ss.str();
}/*}}}*/
template<typename Test> void BenchmarkRunner::executeTest()/*{{{*/
{
    m_threadPool.setTestFunction(&Test::run);
    const size_t memorySizeT = m_memorySize * (GiB / sizeof(Scalar));
    Memory const memoryEnd = m_memory + memorySizeT;
    for (const size_t size : Test::sizes()) {
        if (size == 0) {
            continue;
        }
        const size_t sizeT = size / sizeof(Scalar);
        std::stringstream ss0;
        ss0 << Test::name() << " (" << prettyBytes(size) << ')';
        const std::string benchName = ss0.str();
        if (m_only.empty() || m_only == benchName) {
            for (Memory m = m_memory; m + Test::sliceSizeT() <= memoryEnd; m += Test::stride()) {
                for (const size_t offsetPerThread : Test::offsetsPerThread()) {
                    Benchmark::setColumnData("offset per thread", toString(offsetPerThread));
                    Benchmark::setColumnData("memory location", toString((m - m_memory) / Test::stride()));
                    const int repetitions = std::max<int>(1, Test::sliceSizeT() / sizeT);
                    Benchmark bench(benchName, size * repetitions * m_threadCount * Test::interpretFactor(), Test::interpretUnit());
                    Timer timer;
                    size_t offset = 0;
                    m_threadPool.executeWith(m, [&offset, size, offsetPerThread] { return offset += offsetPerThread; }, sizeT, repetitions);
                    Test::run({m, &timer, 0, sizeT, repetitions});
                    m_threadPool.waitReady();
                    bench.addTiming(timer);
                    m_threadPool.eachTimer([&bench](const Timer &t) { bench.addTiming(t); });
                    bench.Print();
                    if (m_threadCount == 1) {
                        break; // with just a single thread, the thread offset makes no difference
                    }
                }
            }
        }
    }
}/*}}}*/
BenchmarkRunner::BenchmarkRunner()/*{{{*/
    : m_coreIds(parseOnlyCpus()),
    m_threadCount(1),
    m_maxMemorySize(largestMemorySize() / GiB),
    m_memorySize(valueForArgument("--size", m_maxMemorySize)),
    m_only(valueForArgument("--only", std::string())),
    m_memory(nullptr)
{
    Benchmark::addColumn("memory location");
    Benchmark::addColumn("offset per thread");
#ifndef NO_LIBNUMA/*{{{*/
    if (numa_available() == -1) {
        std::cerr << "NUMA interface does not work. Abort." << std::endl;
        return;
    }
#ifdef LINK_STATICALLY
    numa_init();
#endif
#endif/*}}}*/

    if (m_memorySize < 1) {/*{{{*/
        std::cerr << "Need at least 1GB." << std::endl;
        return;
    }
    if (m_memorySize > m_maxMemorySize) {
        std::cerr << "Not enough memory available. Expect crashes/OOM kills." << std::endl;
    }/*}}}*/
    m_memory = Vc::malloc<Scalar, Vc::AlignOnPage>(m_memorySize * GiB / sizeof(Scalar));
#ifndef NO_LIBNUMA/*: mbind {{{*/
    struct bitmask *nodemask = numa_get_mems_allowed();
    std::cerr << " numa_max_node: " << numa_max_node()
        << "\n numa_max_possible_node: " << numa_max_possible_node()
        << "\n numa_num_configured_nodes: " << numa_num_configured_nodes()
        << "\n numa_num_configured_cpus: " << numa_num_configured_cpus()
        << "\n numa_preferred: " << numa_preferred()
        << "\n numa_pagesize: " << numa_pagesize()
        << "\n numa_get_mems_allowed: " << nodemask->maskp[0]
        << std::endl;
    char *addr = reinterpret_cast<char *>(m_memory);
    char *const endAddr = addr + (m_memorySize * GiB);
    const unsigned long maxnode = numa_max_node() + 1;
    size_t bindStride = std::max(1ul, m_memorySize / maxnode) * GiB;
    unsigned long mbindMask[2] = { 0ul, 0ul };
    for (size_t i = 0; i < maxnode; ++i) {
        if (numa_bitmask_isbitset(nodemask, i)) {
            if (addr + bindStride > endAddr) {
                bindStride = endAddr - addr;
            }
            mbindMask[0] = 1 << i;
            //std::cerr << "calling mbind(" << (void*)addr << ", 0x"
            //          << std::hex << bindStride << std::dec << ", MPOL_BIND, "
            //          << mbindMask << ", " << std::min(2ul, maxnode)
            //          << ", MPOL_MF_STRICT | MPOL_MF_MOVE)\n";
            if (-1 == mbind(addr, bindStride, MPOL_BIND, &mbindMask[0],
                        std::min(2ul, maxnode), MPOL_MF_STRICT | MPOL_MF_MOVE)) {
                std::cerr << "mbind failed for node " << i << " of " <<
                   maxnode - 1 << " with " << strerror(errno) << std::endl;
                return;
            } else {
                std::cout << "pinned " << bindStride << " Bytes of memory to NUMA node " << i << std::endl;
            }
            addr += bindStride;
        }
    }
    /*
       EINVAL An invalid value was specified for flags or mode;
              or addr + len was less than addr;
              or addr is not a multiple of the system page size.
              Or, mode is MPOL_DEFAULT and nodemask specified a nonempty set;
              Or, maxnode exceeds a kernel-imposed limit.
              Or, the mode argument specified both MPOL_F_STATIC_NODES and MPOL_F_RELATIVE_NODES.

              or mode is MPOL_BIND or MPOL_INTERLEAVE and nodemask is empty.
              Or, nodemask specifies one or more node IDs that are greater than the maximum supported node ID.
              Or, none of the node IDs specified by nodemask are on-line and allowed by the process's current cpuset context,
              or none of the specified nodes contain memory.
    */
#endif/*}}}*/
    mlockall(MCL_CURRENT);

    if (!m_coreIds.empty()) {/*{{{*/
        Benchmark::addColumn("CPU_ID");
        for (auto cpuRange : m_coreIds) {
            const int cpuid = cpuRange.first;
            cpuRange.first += cpuRange.step;
            std::ostringstream str;
            str << cpuid;
            m_threadCount = 1;
            for (int id = cpuRange.first; id <= cpuRange.last; id += cpuRange.step) {
                str << ',' << id;
                ++m_threadCount;
            }
            Benchmark::setColumnData("CPU_ID", str.str());
            pinToCpu(cpuid);
            m_threadPool.setPinning(cpuRange);
            executeAllTests();
        }
        return;
    }
/*}}}*/
#ifdef NO_LIBNUMA/*{{{*/
    cpu_set_t cpumask;
    sched_getaffinity(0, sizeof(cpu_set_t), &cpumask);
    int cpucount = cpuCount(&cpumask);
    Benchmark::addColumn("CPU_ID");
    for (int cpuid = valueForArgument("--firstCpu", 1); cpuid < cpucount; cpuid += valueForArgument("--cpuStep", 6)) {
        std::ostringstream str;
        str << cpuid;
        Benchmark::setColumnData("CPU_ID", str.str());
        pinToCpu(cpuid);
        executeAllTests();
    }/*}}}*/
#else/*{{{*/
    // libnuma defines:
    // node: an area where all memory as the same speed as seen from a particular CPU. A node
    //       can contain multiple CPUs
    // cpu: a hardware thread?
    int nodeCount = numa_max_node();
    struct bitmask *nodemask = 0;
    if (nodeCount < 0) {
        std::cerr << "libnuma does not report any NUMA nodes\n";
        nodeCount = 0;
    } else {
        nodemask = numa_allocate_nodemask();
    }
    Benchmark::addColumn("NUMA_ID");
    for (int numaId = valueForArgument("--firstNode", 0); numaId <= nodeCount; numaId += valueForArgument("--nodeStep", 1)) {
        std::ostringstream str;
        str << numaId;
        Benchmark::setColumnData("NUMA_ID", str.str());

        if (nodemask) {
            numa_bitmask_clearall(nodemask);
            numa_bitmask_setbit(nodemask, numaId);
            numa_bind(nodemask);
        }
        executeAllTests();
    }
#endif/*}}}*/
#ifndef NO_LIBNUMA/*{{{*/
    if (nodemask) {
        numa_free_nodemask(nodemask);
    }
#endif/*}}}*/
}/*}}}*/
void BenchmarkRunner::executeAllTests()/*{{{*/
{
    //executeTest<TestBzero>();
    executeTest<TestRead>();
    executeTest<TestReadPrefetch>();
    executeTest<TestAddOneStrided>();
    executeTest<TestAddOne>();
    executeTest<TestAddOnePrefetch>();
    executeTest<TestWrite>();
    executeTest<TestStream>();
    //executeTest<TestReadLatency>();
    Benchmark::finalize();
}
/*}}}*/
template<typename Test> static void printDocumentation()/*{{{*/
{
    std::cout
        << "********************************************************************************\n" << Test::name()
        << "\n\n" << Test::description()
        << "\n\nslice size: " << prettyBytes(Test::sliceSizeT() * sizeof(Scalar))
        << "\n     sizes: ";
    for (const size_t s : Test::sizes()) {
        std::cout << prettyBytes(s) << ' ';
    }
    std::cout << "\n   offsets: ";
    for (const size_t s : Test::offsetsPerThread()) {
        std::cout << prettyBytes(s * sizeof(Scalar)) << ' ';
    }
    std::cout << "\n    stride: " << prettyBytes(Test::stride() * sizeof(Scalar)) << std::endl;
}/*}}}*/
int bmain()/*{{{*/
{
    CpuId::init();
    ArgumentVector::iterator it = std::find(g_arguments.begin(), g_arguments.end(), "--listTests");
    if (it != g_arguments.end()) {
        //printDocumentation<TestBzero>();
        printDocumentation<TestRead>();
        printDocumentation<TestReadPrefetch>();
        printDocumentation<TestAddOneStrided>();
        printDocumentation<TestAddOne>();
        printDocumentation<TestAddOnePrefetch>();
        printDocumentation<TestWrite>();
        printDocumentation<TestStream>();
        //printDocumentation<TestReadLatency>();
        return 0;
    }
    BenchmarkRunner runner;
    return 0;
}/*}}}*/

// vim: foldmethod=marker
