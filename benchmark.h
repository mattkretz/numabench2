/*  This file is part of the Vc library.

    Copyright (C) 2009 Matthias Kretz <kretz@kde.org>

    Vc is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Vc is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Vc.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <list>
#include <vector>
#include <algorithm>
#include <time.h>
#include <cstring>
#include <string>
#include <fstream>
#ifndef VC_BENCHMARK_NO_MLOCK
#include <sys/mman.h>
#endif
#ifdef _MSC_VER
#include <windows.h>
#include <float.h>
#else
#include <cmath>
#endif
#ifdef __APPLE__
#include <mach/mach_time.h>
// method to get monotonic mac time, inspired by 
// http://www.wand.net.nz/~smr26/wordpress/2009/01/19/monotonic-time-in-mac-os-x/
#endif

#include "tsc.h"

#ifndef _MSC_VER
static inline double convertTimeSpec(const struct timespec &ts)
{
    return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) * 1e-9;
}
#endif

class Timer
{
private:
#ifdef _MSC_VER
    __int64 m_startTime;
#elif defined(__APPLE__)
    uint64_t m_startTime;
#else
    struct timespec m_startTime;
#endif
    TimeStampCounter m_tsc;

    struct ElapsedTime
    {
        double real;
        double cycles;
    } m_elapsed;

public:
    inline void start() {
#ifdef _MSC_VER
        QueryPerformanceCounter((LARGE_INTEGER *)&m_startTime);
#elif defined(__APPLE__)
        m_startTime = mach_absolute_time();
#else
        clock_gettime( CLOCK_MONOTONIC, &m_startTime );
#endif
        m_tsc.start();
    }
    inline void stop()
    {
        m_tsc.stop();
#ifdef _MSC_VER
        __int64 real = 0, freq = 0;
        QueryPerformanceCounter((LARGE_INTEGER *)&real);
        QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
        m_elapsed.real = static_cast<double>(real - m_startTime) / freq;
#elif defined(__APPLE__)
        uint64_t real = mach_absolute_time();
        static mach_timebase_info_data_t info = {0,0};

        if (info.denom == 0)
            mach_timebase_info(&info);

        uint64_t nanos = (real - m_startTime ) * (info.numer / info.denom);
        m_elapsed.real = nanos * 1e-9;
#else
        struct timespec real;
        clock_gettime( CLOCK_MONOTONIC, &real );
        m_elapsed.real = convertTimeSpec(real) - convertTimeSpec(m_startTime);
#endif
        m_elapsed.cycles = m_tsc.cycles();
    }
    inline double realTime() const {
        return m_elapsed.real;
    }
    inline double cycles() const {
        return m_elapsed.cycles;
    }
};

class Benchmark
{
    friend int main(int, char**);
    class FileWriter
    {
        public:
            FileWriter(const std::string &filename);
            ~FileWriter();

            void declareData(const std::string &name, const std::list<std::string> &header);
            void addDataLine(const std::list<std::string> &data);

            void addColumn(const std::string &name);
            void setColumnData(const std::string &name, const std::string &data);
            void finalize() { m_finalized = true; }

        private:
            std::ofstream m_file;
            bool m_finalized;
            std::string m_currentName;
            std::list<std::string> m_header;
            struct ExtraColumn
            {
                ExtraColumn(const std::string &n) : name(n) {}
                std::string name;
                std::string data;
                inline bool operator==(const ExtraColumn &rhs) const { return name == rhs.name; }
                inline bool operator==(const std::string &rhs) const { return name == rhs; }
            };
            std::list<ExtraColumn> m_extraColumns;
    };
public:
    static inline void addColumn(const std::string &name) { if (s_fileWriter) s_fileWriter->addColumn(name); }
    static inline void setColumnData(const std::string &name, const std::string &data) {
        if (s_fileWriter) {
            s_fileWriter->setColumnData(name, data);
        } else {
            std::cout << "Benchmarking " << name << " " << data << std::endl;
        }
    }
    static inline void finalize() { if (s_fileWriter) s_fileWriter->finalize(); }

    explicit Benchmark(const std::string &name, double factor = 0., const std::string &X = std::string());
    void changeInterpretation(double factor, const char *X)
    {
        fFactor = factor;
        fX = X;
    }

    void addTiming(const Timer &t);
    void Print();

private:
    void printMiddleLine() const;
    void printBottomLine() const;

    const std::string fName;
    double fFactor;
    std::string fX;
    double m_mean[2];
    double m_stddev[2];
    int m_dataPointsCount;
    static FileWriter *s_fileWriter;
};

Benchmark::FileWriter::FileWriter(const std::string &filename)
    : m_finalized(false)
{
    std::string fn = filename;
    m_file.open(fn.c_str());

    if (Benchmark::s_fileWriter == 0) {
        Benchmark::s_fileWriter = this;
    }
}

Benchmark::FileWriter::~FileWriter()
{
    if (m_file.is_open()) {
        m_file.flush();
        m_file.close();
    }
    if (Benchmark::s_fileWriter == this) {
        Benchmark::s_fileWriter = 0;
    }
}

void Benchmark::FileWriter::declareData(const std::string &name, const std::list<std::string> &header)
{
    m_currentName = '"' + name + '"';
    if (m_header != header) {
        if (m_header.empty()) {
            m_file << "Version 3\n";
        }
        m_header = header;
        m_file << "\"benchmark.name\"\t\"benchmark.arch\"";
        for (std::list<ExtraColumn>::const_iterator i = m_extraColumns.begin();
                i != m_extraColumns.end(); ++i) {
            m_file << "\t\"" << i->name << '"';
        }
        for (std::list<std::string>::const_iterator i = header.begin();
                i != header.end(); ++i) {
            m_file << '\t' << *i;
        }
        m_file << "\n";
    }
}

void Benchmark::FileWriter::addDataLine(const std::list<std::string> &data)
{
    m_file << m_currentName << "\t\"";
    switch (static_cast<Vc::Implementation>(Vc::CurrentImplementation::Implementation)) {
    case Vc::ScalarImpl: m_file << "Scalar"; break;
    case Vc::SSE2Impl:   m_file << "SSE2";   break;
    case Vc::SSE3Impl:   m_file << "SSE3";   break;
    case Vc::SSSE3Impl:  m_file << "SSSE3";  break;
    case Vc::SSE41Impl:  m_file << "SSE4.1"; break;
    case Vc::SSE42Impl:  m_file << "SSE4.2"; break;
    case Vc::AVXImpl:    m_file << "AVX";    break;
    case Vc::AVX2Impl:    m_file << "AVX2";  break;
    }

    const auto extraInstructions = static_cast<Vc::ExtraInstructions>(Vc::CurrentImplementation::ExtraInstructions);
    if (extraInstructions & Vc::Sse4aInstructions ) m_file << "+SSE4a";
    if (extraInstructions & Vc::XopInstructions   ) m_file << "+XOP";
    if (extraInstructions & Vc::Fma4Instructions  ) m_file << "+FMA4";
    if (extraInstructions & Vc::PopcntInstructions) m_file << "+POPCNT";
    if (extraInstructions & Vc::FmaInstructions   ) m_file << "+FMA";
    m_file << '"';

    for (std::list<ExtraColumn>::const_iterator i = m_extraColumns.begin();
            i != m_extraColumns.end(); ++i) {
        m_file << '\t' << i->data;
    }
    for (std::list<std::string>::const_iterator i = data.begin();
            i != data.end(); ++i) {
        m_file << '\t' << *i;
    }
    m_file << "\n";
}

void Benchmark::FileWriter::addColumn(const std::string &name)
{
    if (!m_finalized) {
        if (m_header.empty()) {
            if (m_extraColumns.end() == std::find(m_extraColumns.begin(), m_extraColumns.end(), name)) {
                m_extraColumns.push_back(name);
            }
        } else {
            std::cerr << "call addColumn before the first benchmark prints its data" << std::endl;
        }
    }
}

void Benchmark::FileWriter::setColumnData(const std::string &name, const std::string &data)
{
    for (std::list<ExtraColumn>::iterator i = m_extraColumns.begin();
            i != m_extraColumns.end(); ++i) {
        if (*i == name) {
            i->data = '"' + data + '"';
            break;
        }
    }
}

Benchmark::FileWriter *Benchmark::s_fileWriter = 0;

Benchmark::Benchmark(const std::string &_name, double factor, const std::string &X)
    : fName(_name), fFactor(factor), fX(X), m_dataPointsCount(0)
{
    for (int i = 0; i < 2; ++i) {
        m_mean[i] = m_stddev[i] = 0.;
    }
    enum {
        WCHARSIZE = sizeof("━") - 1
    };
    if (!s_fileWriter) {
        const bool interpret = (fFactor != 0.);
        char header[128 * WCHARSIZE];
        std::memset(header, 0, 128 * WCHARSIZE);
        std::strcpy(header,
                "┏━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┓");
        if (!interpret) {
            header[(69 - 17) * WCHARSIZE] = '\0';
        }
        const int titleLen = fName.length();
        const int headerLen = std::strlen(header) / WCHARSIZE;
        int offset = (headerLen - titleLen) / 2;
        if (offset > 0) {
            --offset;
            std::string name = ' ' + fName + ' ';
            char *ptr = &header[offset * WCHARSIZE];
            std::memcpy(ptr, name.c_str(), name.length());
            std::memmove(ptr + name.length(), ptr + name.length() * WCHARSIZE, (headerLen - offset - name.length()) * WCHARSIZE + 1);
            std::cout << header << std::flush;
        } else {
            std::cout << fName << std::flush;
        }
    }
}

static inline void prettyPrintSeconds(double v)
{
    static const char prefix[] = { ' ', 'm', 'u', 'n', 'p' };
    if (v == 0.) {
        std::cout << "      0       ";
    } else if (v < 2.) {
        int i = 0;
        do {
            v *= 1000.;
            ++i;
        } while (v < 1.);
        std::cout << std::setw(11) << v << ' ' << prefix[i] << 's';
    } else if (v > 60.) {
        std::cout << std::setw(10) << v / 60. << " min";
    } else {
        std::cout << std::setw(12) << v << " s";
    }
}

#ifdef isfinite
#undef isfinite
#endif

static inline void prettyPrintCount(double v)
{
    static const char prefix[] = { ' ', 'k', 'M', 'G', 'T', 'P', 'E' };
    int i = 0;
#ifdef _MSC_VER
    if (_finite(v)) {
#elif defined(__INTEL_COMPILER)
    if (::isfinite(v)) {
#else
    if (std::isfinite(v)) {
#endif
        if (v < 1000.) {
            std::cout << std::setw(14) << v;
            return;
        }
        while (v >= 1000.) {
            v *= 0.001;
            ++i;
        }
    }
    std::cout << std::setw(12) << v << ' ' << prefix[i];
}

static inline void prettyPrintError(double v)
{
    std::stringstream ss;
    ss << "± " << v << " %";
    std::cout << std::setw(15) << ss.str();
}

static inline std::list<std::string> &operator<<(std::list<std::string> &list, const std::string &data)
{
    std::ostringstream str;
    str << '"' << data << '"';
    list.push_back(str.str());
    return list;
}

static inline std::list<std::string> &operator<<(std::list<std::string> &list, const char *data)
{
    std::ostringstream str;
    str << '"' << data << '"';
    list.push_back(str.str());
    return list;
}

static inline std::list<std::string> &operator<<(std::list<std::string> &list, double data)
{
    std::ostringstream str;
    str << data;
    list.push_back(str.str());
    return list;
}

static inline std::list<std::string> &operator<<(std::list<std::string> &list, int data)
{
    std::ostringstream str;
    str << data;
    list.push_back(str.str());
    return list;
}

static inline std::list<std::string> &operator<<(std::list<std::string> &list, unsigned int data)
{
    std::ostringstream str;
    str << data;
    list.push_back(str.str());
    return list;
}

static inline std::list<std::string> &operator<<(std::list<std::string> &list, unsigned long long data)
{
    std::ostringstream str;
    str << data;
    list.push_back(str.str());
    return list;
}

static std::string centered(const std::string &s, const int size = 16)
{
    const int missing = size - s.length();
    if (missing < 0) {
        return s.substr(0, size);
    } else if (missing == 0) {
        return s;
    }
    const int left = missing - missing / 2;
    std::string r(size, ' ');
    r.replace(left, s.length(), s);
    return r;
}

inline void Benchmark::printMiddleLine() const
{
    const bool interpret = (fFactor != 0.);
    std::cout << "\n"
        "┠────────────────╂────────────────"
        << (interpret ?
                "╂────────────────╂────────────────╂────────────────┨" : "┨");
}
inline void Benchmark::printBottomLine() const
{
    const bool interpret = (fFactor != 0.);
    std::cout << "\n"
        "┗━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━"
        << (interpret ?
                "┻━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━┛" : "┛") << std::endl;
}

inline void Benchmark::addTiming(const Timer &t) {
    m_mean[0] += t.realTime();
    m_mean[1] += t.cycles();
    m_stddev[0] += t.realTime() * t.realTime();
    m_stddev[1] += t.cycles() * t.cycles();
    ++m_dataPointsCount;
}

inline void Benchmark::Print()
{
    std::streambuf *backup = std::cout.rdbuf();
    if (s_fileWriter) {
        std::cout.rdbuf(0);
    }
    const bool interpret = (fFactor != 0.);

    std::list<std::string> header;
    header
        << "Real_time" << "Real_time_stddev"
        << "Cycles" << "Cycles_stddev"
    ;

    // ┃ ━ ┏ ┓ ┗ ┛ ┣ ┫ ┳ ┻ ╋ ┠ ─ ╂ ┨
    std::cout << "\n"
        << "┃   Real time    ┃     Cycles     ┃";
    if (interpret) {
        std::cout << centered(fX + "/s [Real]") << "┃";
        std::cout << centered(fX + "/cycle")    << "┃";
        std::cout << centered("cycles/" + fX)   << "┃";
        std::string X = fX;
        for (unsigned int i = 0; i < X.length(); ++i) {
            if (X[i] == ' ') {
                X[i] = '_';
            }
        }
        header
            << X + "/Real_time" << X + "/Real_time_stddev"
            << X + "/Cycles" << X + "/Cycles_stddev"
            << "number_of_" + X;
    }
    printMiddleLine();
    if (s_fileWriter) {
        s_fileWriter->declareData(fName, header);
    }

    const double normalization = 1. / m_dataPointsCount;

    std::list<std::string> dataLine;
    m_mean[0] *= normalization;
    m_stddev[0] = std::sqrt(m_stddev[0] * normalization - m_mean[0] * m_mean[0]);
    dataLine << m_mean[0] << m_stddev[0];
    m_mean[1] *= normalization;
    m_stddev[1] = std::sqrt(m_stddev[1] * normalization - m_mean[1] * m_mean[1]);
    dataLine << m_mean[1] << m_stddev[1];
    double stddevint[2];
    stddevint[0] = fFactor * m_stddev[0] / (m_mean[0] * m_mean[0]);
    stddevint[1] = fFactor * m_stddev[1] / (m_mean[1] * m_mean[1]);
    dataLine << fFactor / m_mean[0] << stddevint[0];
    dataLine << fFactor / m_mean[1] << stddevint[1];
    dataLine << fFactor;

    std::cout << "\n┃ ";
    prettyPrintSeconds(m_mean[0]);
    std::cout << " ┃ ";
    prettyPrintCount(m_mean[1]);
    std::cout << " ┃ ";
    if (interpret) {
        prettyPrintCount(fFactor / m_mean[0]);
        std::cout << " ┃ ";
        prettyPrintCount(fFactor / m_mean[1]);
        std::cout << " ┃ ";
        prettyPrintCount(m_mean[1] / fFactor);
        std::cout << " ┃ ";
    }
    std::cout << "\n┃ ";
    prettyPrintError(m_stddev[0] * 100. / m_mean[0]);
    std::cout << " ┃ ";
    prettyPrintError(m_stddev[1] * 100. / m_mean[1]);
    std::cout << " ┃ ";
    if (interpret) {
        prettyPrintError(m_stddev[0] * 100. / m_mean[0]);
        std::cout << " ┃ ";
        prettyPrintError(m_stddev[1] * 100. / m_mean[1]);
        std::cout << " ┃ ";
        prettyPrintError(m_stddev[1] * 100. / m_mean[1]);
        std::cout << " ┃ ";
    }
    printBottomLine();
    if (s_fileWriter) {
        s_fileWriter->addDataLine(dataLine);
        std::cout.rdbuf(backup);
    }
}

typedef std::vector<std::string> ArgumentVector;
ArgumentVector g_arguments;

int bmain();
const char *printHelp2 =
"  -t <seconds>        maximum time to run a single benchmark (10s)\n"
"  -cpu (all|any|<id>) CPU to pin the benchmark to\n"
"                      all: test every CPU id in sequence\n"
"                      any: don't pin and let the OS schedule\n"
"                      <id>: pin to the specific CPU\n";
#define SET_HELP_TEXT(str) \
    int _set_help_text_init() { \
        printHelp2 = str; \
        return 0; \
    } \
    int _set_help_text_init_ = _set_help_text_init()

#include "cpuset.h"

void printHelp(const char *name) {
    std::cout << "Usage " << name << " [OPTION]...\n"
        << "Measure throughput and latency of memory in steps of 1GB\n\n"
        << "  -h, --help          print this message\n"
        << "  -o <filename>       output measurements to a file instead of stdout\n";
    if (printHelp2) {
        std::cout << printHelp2;
    }
    std::cout << "\nReport bugs to vc-devel@compeng.uni-frankfurt.de\n"
        << "Vc Homepage: http://compeng.uni-frankfurt.de/index.php?id=Vc\n"
        << std::flush;
}

int main(int argc, char **argv)
{
#ifdef SCHED_FIFO_BENCHMARKS
    if (SCHED_FIFO != sched_getscheduler(0)) {
        // not realtime priority, check whether the rtwrapper executable exists
        execv("./rtwrapper", argv);
        // if the execv call works, great. If it doesn't we just continue, but without realtime prio
    }
#endif

    int i = 2;
    Benchmark::FileWriter *file = 0;
    while (argc > i) {
        if (std::strcmp(argv[i - 1], "-o") == 0) {
            file = new Benchmark::FileWriter(argv[i]);
            i += 2;
        } else if (std::strcmp(argv[i - 1], "--help") == 0 ||
                    std::strcmp(argv[i - 1], "-help") == 0 ||
                    std::strcmp(argv[i - 1], "-h") == 0) {
            printHelp(argv[0]);
            return 0;
        } else {
            g_arguments.push_back(argv[i - 1]);
            ++i;
        }
    }
    if (argc == i) {
        if (std::strcmp(argv[i - 1], "--help") == 0 ||
                std::strcmp(argv[i - 1], "-help") == 0 ||
                std::strcmp(argv[i - 1], "-h") == 0) {
            printHelp(argv[0]);
            return 0;
        }
        g_arguments.push_back(argv[i - 1]);
    }

    int r = bmain();
    Benchmark::finalize();
    delete file;
    return r;
}

#endif // BENCHMARK_H
