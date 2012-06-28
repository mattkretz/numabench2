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
#define NO_LIBNUMA 1
#ifdef NO_LIBNUMA
#include "cpuset.h"
#else
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
struct CpuRange/*{{{*/
{
    int first;
    int last;
    int step;
};/*}}}*/
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

    private:
        void mainLoop() // thread
        {
            m_mutex.lock();
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
                    cpu_set_t cpumask;
                    cpuZero(&cpumask);
                    cpuSet(m_cpuId, &cpumask);
                    sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
                    m_cpuId = -1;
                    continue;
                }


                // do the work
                m_testFunction(m_arguments);
            } while (!m_exit);
            m_waitForEnd.oneReady();
            m_mutex.unlock();
        }
};/*}}}*/
class ThreadPool/*{{{*/
{
    OneWaitsForN m_waitForEnd;
    std::condition_variable_any m_waitForStart;
    std::vector<std::shared_ptr<ThreadData>> m_workers;
    static int maxThreadCount() { cpu_set_t cpumask; sched_getaffinity(0, sizeof(cpu_set_t), &cpumask); return cpuCount(&cpumask); }

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
struct TestDefaults/*{{{*/
{
    /// in Bytes
    static std::vector<size_t> sizes() { return { 1 * GiB, CpuId::L3Data() / 2, CpuId::L2Data() / 2, CpuId::L1Data() / 2 }; }
    /// in #Scalars, wholeSize in Bytes
    static constexpr size_t offsetPerThread(size_t wholeSize) {
        return 0;//wholeSize / (sizeof(Scalar) * 60); //< GiB ? 7 * ScalarsInCacheLine : 9 * ScalarsInPage;
    }
    static constexpr double interpretFactor() { return 1.; }
    static constexpr const char *interpretUnit() { return "Byte"; }
    /// in #Scalars
    static constexpr size_t stride() { return GiB / sizeof(Scalar); }
};/*}}}*/
struct TestBzero : public TestDefaults/*{{{*/
{
    static constexpr const char *name() { return "bzero"; }
    static void run(const TestArguments &args)
    {
        args.timer->start();
        for (int rep = 0; rep < args.repetitions; ++rep) {
            bzero(args.mem + args.offset, (args.size - args.offset) * sizeof(Scalar));
            if (args.offset > 0) {
                bzero(args.mem, args.offset * sizeof(Scalar));
            }
        }
        args.timer->stop();
    }
};/*}}}*/
struct TestAddOneStrided : public TestDefaults/*{{{*/
{
    static std::vector<size_t> sizes() { return { 8 * GiB }; }
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
};/*}}}*/
struct TestAddOneStrided2 : public TestDefaults/*{{{*/
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
template<bool Prefetch> struct TestAddOneBase : public TestDefaults/*{{{*/
{
    static void run(const TestArguments &args)
    {
        const Vector one = Vector::One();

        const Memory mStart[2] = { args.mem + args.offset, args.mem  };
        const Memory mEnd  [2] = { args.mem + args.size  , mStart[0] };

        args.timer->start();
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (int i = 0; i < 2; ++i) {
                for (Memory m = mStart[i] + 3 * Vector::Size; m < mEnd[i]; m += 4 * Vector::Size) {
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
template<bool Prefetch> struct TestReadBase : public TestDefaults/*{{{*/
{
    static void run(const TestArguments &args)
    {
        const Memory mStart[2] = { args.mem + args.offset, args.mem  };
        const Memory mEnd  [2] = { args.mem + args.size  , mStart[0] };
        args.timer->start();
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (int i = 0; i < 2; ++i) {
                for (Memory m = mStart[i] + 3 * Vector::Size; m < mEnd[i]; m += 4 * Vector::Size) {
                    if (Prefetch) {
                        Vc::prefetchForOneRead(m + 4032 / sizeof(Scalar));
                    }
                    const Vector v0(m - 3 * Vector::Size);
                    const Vector v1(m - 2 * Vector::Size);
                    const Vector v2(m - 1 * Vector::Size);
                    const Vector v3(m - 0 * Vector::Size);
                    asm("" :: "x"(v0.data()), "x"(v1.data()), "x"(v2.data()), "x"(v3.data()));
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
struct TestReadLatency : public TestDefaults
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
/*SET_HELP_TEXT{{{*/
#ifdef NO_LIBNUMA
SET_HELP_TEXT(
        "  --firstCpu <id>\n"
        "  --cpuStep <id>\n"
        "  --size <GiB>\n"
        "  --only <test function>\n"
        "  --cores <firstId-lastId[:step][,firstId-lastId[:step][...]]>\n"
        );
#else
SET_HELP_TEXT(
        "  --firstNode <id>\n"
        "  --nodeStep <id>\n"
        "  --size <GiB>\n"
        "  --only <test function>\n"
        "  --cores <firstId-lastId[:step][,firstId-lastId[:step][...]]>\n"
        );
#endif/*}}}*/
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
        if (m_only.empty() || m_only == ss0.str()) {
            for (Memory m = m_memory; m + sizeT <= memoryEnd; m += Test::stride()) {
                std::stringstream ss;
                ss << ss0.str();
                ss << " @ " << prettyBytes((m - m_memory) * sizeof(Scalar));
                const int repetitions = std::max<int>(1, Test::stride() / sizeT);
                Benchmark bench(ss.str().c_str(), size * repetitions * m_threadCount * Test::interpretFactor(), Test::interpretUnit());
                Timer timer;
                size_t offset = 0;
                m_threadPool.executeWith(m, [&offset, size] { return offset += Test::offsetPerThread(size); }, sizeT, repetitions);
                Test::run({m, &timer, 0, sizeT, repetitions});
                m_threadPool.waitReady();
                bench.addTiming(timer);
                m_threadPool.eachTimer([&bench](const Timer &t) { bench.addTiming(t); });
                bench.Print();
            }
        }
    }
}/*}}}*/
inline std::ostream &operator<<(std::ostream &out, const CpuRange &range)/*{{{*/
{
    return out << "Range " << range.first << " - " << range.last << ", step " << range.step << '\n';
}
/*}}}*/
std::vector<CpuRange> parseOnlyCpus()/*{{{*/
{
    enum State {
        ReadFirst, ReadLast, ReadStep
    };
    const std::string cpusStrings = valueForArgument("--cores", std::string());
    std::vector<CpuRange> r;
    CpuRange range = { 0, 0, 1 };
    State state = ReadFirst;
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
            state = ReadFirst;
            r.push_back(range);
            range = { 0, 0, 1 };
        }
    }
    r.push_back(range);
    return r;
}
/*}}}*/
BenchmarkRunner::BenchmarkRunner()/*{{{*/
    : m_coreIds(parseOnlyCpus()),
    m_threadCount(1),
    m_maxMemorySize(largestMemorySize() / GiB),
    m_memorySize(valueForArgument("--size", m_maxMemorySize)),
    m_only(valueForArgument("--only", std::string())),
    m_memory(nullptr)
{
#ifndef NO_LIBNUMA/*{{{*/
    if (numa_available() == -1) {
        std::cerr << "NUMA interface does not work. Abort." << std::endl;
        return;
    }
#ifdef LINK_STATICALLY
    numa_init();
#endif

    // first make sure we don't get interleaved memory; this would defeat the purpose of this
    // benchmark
    //numa_set_interleave_mask(numa_no_nodes);
#endif/*}}}*/

    if (m_memorySize < 1) {/*{{{*/
        std::cerr << "Need at least 1GB." << std::endl;
        return;
    }
    if (m_memorySize > m_maxMemorySize) {
        std::cerr << "Not enough memory available. Expect crashes/OOM kills." << std::endl;
    }/*}}}*/
    m_memory = Vc::malloc<Scalar, Vc::AlignOnPage>(m_memorySize * GiB / sizeof(Scalar));
    mlockall(MCL_CURRENT);

    if (!m_coreIds.empty()) {/*{{{*/
        Benchmark::addColumn("CPU_ID");
        for (auto cpuRange : m_coreIds) {
            cpu_set_t cpumask;
            sched_getaffinity(0, sizeof(cpu_set_t), &cpumask);
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
            cpuZero(&cpumask);
            cpuSet(cpuid, &cpumask);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
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
        cpuZero(&cpumask);
        cpuSet(cpuid, &cpumask);
        sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
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
    executeTest<TestBzero>();
    executeTest<TestRead>();
    executeTest<TestReadPrefetch>();
    executeTest<TestAddOneStrided>();
    executeTest<TestAddOne>();
    executeTest<TestAddOnePrefetch>();
    executeTest<TestReadLatency>();
    Benchmark::finalize();
}
/*}}}*/
int bmain()/*{{{*/
{
    CpuId::init();
    BenchmarkRunner runner;
    return 0;
}/*}}}*/

// vim: foldmethod=marker
