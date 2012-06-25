/*
    Copyright (C) 2009 Matthias Kretz <kretz@kde.org>

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) version 3.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public License
    along with this library; see the file COPYING.LIB.  If not, write to
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
    Boston, MA 02110-1301, USA.

*/

#ifndef TSC_H
#define TSC_H

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__rdtsc)
#endif

class TimeStampCounter
{
    public:
        void start();
        void stop();
        unsigned long long cycles() const;

    private:
        union Data {
            unsigned long long a;
            unsigned int b[2];
        } m_start, m_end;
};

inline void TimeStampCounter::start()
{
#ifdef _MSC_VER
    m_start.a = __rdtscp();
#else
    asm volatile("rdtscp" : "=a"(m_start.b[0]), "=d"(m_start.b[1]) :: "ecx" );
#endif
}

inline void TimeStampCounter::stop()
{
#ifdef _MSC_VER
    m_end.a = __rdtscp();
#else
    asm volatile("rdtscp" : "=a"(m_end.b[0]), "=d"(m_end.b[1]) :: "ecx" );
#endif
}

inline unsigned long long TimeStampCounter::cycles() const
{
    return m_end.a - m_start.a;
}

#endif // TSC_H
