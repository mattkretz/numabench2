#include <Vc/Vc>/*{{{*/
#include "benchmark.h"
#include <sys/mman.h>
#include <fstream>
#include <thread>
#include <mutex>
#include <atomic>
#include <Vc/cpuid.h>
#include "cpuset.h"

using namespace Vc;/*}}}*/

using std::mutex;
using std::unique_lock;

static int maxThreadCount()
{
    cpu_set_t cpumask;
    sched_getaffinity(0, sizeof(cpu_set_t), &cpumask);
    return cpuCount(&cpumask);
}

class Add1Sweep
{
public:
  Add1Sweep();
};

Add1Sweep::Add1Sweep()
{
  float *__restrict__ mem = Vc::malloc<float, Vc::AlignOnPage>(1024 * 1024 * 1024);
  mlockall(MCL_CURRENT);
  const int maxThreadCount = ::maxThreadCount();

  const size_t counts[] = {
    CpuId::L1Data() / sizeof(float),
    CpuId::L2Data() / sizeof(float),
    CpuId::L3Data() / sizeof(float),
    1024 * 1024 * 1024 / sizeof(float)
  };

  TimeStampCounter tsc;
  const float_v one(1.f);

  for (size_t count : counts) {
    const int iterations = *(std::end(counts) - 1) / count;
    for (int cpu = 0; cpu < maxThreadCount; ++cpu) {
      pinToCpu(cpu);
      for (size_t strideOffset : { 64 / sizeof(float), size_t(0) }) {
        for (size_t stride0 = float_v::Size; stride0 <= 256 * 1024 * 1024; stride0 *= 2) {
          size_t stride = stride0 + strideOffset;
          if (stride >= count) {
            break;
          }
          tsc.start();
          for (int it = 0; it < iterations; ++it) {
            for (size_t start = 0; start < stride; start += float_v::Size) {
              size_t i = start;
              for (; i + 3 * stride < count; i += 4 * stride) {
                const float_v tmp0 = one + float_v(&mem[i + 0 * stride]);
                const float_v tmp1 = one + float_v(&mem[i + 1 * stride]);
                const float_v tmp2 = one + float_v(&mem[i + 2 * stride]);
                const float_v tmp3 = one + float_v(&mem[i + 3 * stride]);
                tmp0.store(&mem[i + 0 * stride]);
                tmp1.store(&mem[i + 1 * stride]);
                tmp2.store(&mem[i + 2 * stride]);
                tmp3.store(&mem[i + 3 * stride]);
              }
              for (; i < count; i += stride) {
                const float_v tmp0 = one + float_v(&mem[i]);
                tmp0.store(&mem[i]);
              }
            }
          }
          tsc.stop();
          double cycles = tsc.cycles();
          cycles /= iterations * count;

          std::cout << cpu << '\t' << count * sizeof(float) / 1024 << '\t' << cycles << '\t' << stride << '\n';
        }
      }
    }
  }
}

class CachePingPong
{
    std::mutex m_mutex;
    std::condition_variable m_wait;
    int m_state;
    double m_cycles;
    double *__restrict__ mem;
    int m_iterations;
    static int maxThreadCount() { cpu_set_t cpumask; sched_getaffinity(0, sizeof(cpu_set_t), &cpumask); return cpuCount(&cpumask); }

public:
    CachePingPong()
        : m_state(0),
        mem(new double[1024 * 1024 * 64])
    {
        mlockall(MCL_CURRENT);

        const size_t counts[] = {
            //CpuId::L1Data() / sizeof(double),
            //CpuId::L2Data() / sizeof(double),
            2 * CpuId::L3Data() / sizeof(double)
        };
        for (size_t count : counts) {
            m_iterations = 16000 * 1024 / count;
            const int benchCpu = 0;
            for (int poisonCpu = 0; poisonCpu < maxThreadCount(); ++poisonCpu) {
                m_state = 0;
                m_cycles = 0.;
                std::thread poison(poisonThreadCaller, this, poisonCpu, count);
                std::thread bench(benchmarkThreadCaller, this, benchCpu, count);
                poison.join();
                bench.join();
                std::cout << poisonCpu << /*'\t' << benchCpu <<*/ '\t' << count * sizeof(double) / 1024
                    << '\t' << m_cycles / m_iterations << '\n';
            }
        }
    }

    static void poisonThreadCaller(CachePingPong *that, int cpuId, size_t count)
    {
        that->poisonThread(cpuId, count);
    }
    static void benchmarkThreadCaller(CachePingPong *that, int cpuId, size_t count)
    {
        that->benchmarkThread(cpuId, count);
    }

    void poisonThread(int cpuId, size_t count)
    {
        pinToCpu(cpuId);
        while (m_state < 2 * m_iterations) {
            unique_lock<mutex> lock(m_mutex);
            if ((m_state & 1) == 1) {
                m_wait.wait(lock);
            }
            // mark the cachelines as modified
            for (size_t i = 0; i < count; i += 4) {
                mem[i] += 1.;
            }
            ++m_state;
            m_wait.notify_one();
        }
    }

    void benchmarkThread(int cpuId, size_t count)
    {
        pinToCpu(cpuId);
        TimeStampCounter tsc;
        const double_v two(2.);
        while (m_state < 2 * m_iterations) {
            unique_lock<mutex> lock(m_mutex);
            if ((m_state & 1) == 0) {
                m_wait.wait(lock);
            }
            tsc.start();
            for (size_t i = 0; i < count; i += 16) {
                const double_v tmp0 = two + double_v(&mem[i +  0]);
                const double_v tmp2 = two + double_v(&mem[i +  4]);
                const double_v tmp4 = two + double_v(&mem[i +  8]);
                const double_v tmp6 = two + double_v(&mem[i + 12]);
                tmp0.store(&mem[i +  0]);
                tmp2.store(&mem[i +  4]);
                tmp4.store(&mem[i +  8]);
                tmp6.store(&mem[i + 12]);
            }
            tsc.stop();
            ++m_state;
            m_wait.notify_one();
            m_cycles += tsc.cycles();
        }
    }
};

int bmain()/*{{{*/
{
    CpuId::init();
    //CachePingPong();
    Add1Sweep();
    return 0;
}
