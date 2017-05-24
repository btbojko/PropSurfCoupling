//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_VERSION_H
#define ANTIOCH_VERSION_H

// C++
#include <iostream>

#define ANTIOCH_MAJOR_VERSION  0
#define ANTIOCH_MINOR_VERSION  5
#define ANTIOCH_MICRO_VERSION  0

#define ANTIOCH_VERSION_AT_LEAST(major, minor, micro) (ANTIOCH_MAJOR_VERSION > (major) || (ANTIOCH_MAJOR_VERSION == (major) && (ANTIOCH_MINOR_VERSION > (minor) || (ANTIOCH_MINOR_VERSION == (minor) && ANTIOCH_MICRO_VERSION > (micro)))))

#define ANTIOCH_BUILD_USER     "bbojko"
#define ANTIOCH_BUILD_ARCH     "x86_64-unknown-linux-gnu"
#define ANTIOCH_BUILD_HOST     "first-rule.wd.navair-rdte.navy.mil"
#define ANTIOCH_BUILD_DATE     "@BUILD_DATE@"
#define ANTIOCH_BUILD_VERSION  "a8a1f9ffe54a51b3703e5b20190d2d361615d5d7"

#define ANTIOCH_LIB_VERSION    "0.5.0"
#define ANTIOCH_LIB_RELEASE    "Development Build"

#define ANTIOCH_CXX            "g++"
#define ANTIOCH_CXXFLAGS       "-g -O2  -Wall -std=c++11"

namespace Antioch
{
  int get_antioch_version();

  template< typename CharT, typename Traits >
  std::basic_ostream<CharT,Traits>&
  antioch_version(std::basic_ostream<CharT,Traits> &os)
  {
    // A little automatic C-style string concatenation goes a long way.
    // It also lets using strings(1) on a binary show something useful.
    return
    os << "--------------------------------------------------------\n"
          "antioch Package: Version = " ANTIOCH_LIB_VERSION " ("
       << get_antioch_version() << ")\n"
          "\n"
          ANTIOCH_LIB_RELEASE "\n"
          "\n"
          "Build Date   = " ANTIOCH_BUILD_DATE     "\n"
          "Build Host   = " ANTIOCH_BUILD_HOST     "\n"
          "Build User   = " ANTIOCH_BUILD_USER     "\n"
          "Build Arch   = " ANTIOCH_BUILD_ARCH     "\n"
          "Build Rev    = " ANTIOCH_BUILD_VERSION  "\n"
          "\n"
          "C++ Config   = " ANTIOCH_CXX " " ANTIOCH_CXXFLAGS "\n"
          "--------------------------------------------------------\n";
  }

  inline void antioch_version_stdout()
  {
    antioch_version(std::cout).flush();
  }

} // end namespace Antioch

#endif //ANTIOCH_VERISON_H
