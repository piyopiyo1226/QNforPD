#pragma once
#ifdef _MSC_VER
#pragma once
#endif

#ifndef ARCSIM_API_HPP
#define ARCSIM_API_HPP

// -----------------------------------------------------------------------------
// API export macro
// -----------------------------------------------------------------------------

#if (defined(WIN32) || defined(_WIN32) || defined(WINCE) || defined(__CYGWIN__))
#   if defined(ARCSIM_API_EXPORT)
#       define ARCSIM_EXPORTS __declspec(dllexport)
#       define ARCSIM_IMPORTS
#   else
#       define ARCSIM_EXPORTS
#       define ARCSIM_IMPORTS __declspec(dllimport)
#   endif
#   define MEM_ALIGN(n) alignas(n)
#elif defined(__GNUC__) && __GNUC__ >= 4
#   define ARCSIM_EXPORTS __attribute__((visibility ("default")))
#   define ARCSIM_IMPORTS __attribute__((visibility ("hidden")))
#   define MEM_ALIGN(n) __attribute__((aligned(n)))
#else
#   define ARCSIM_EXPORTS
#   define ARCSIM_IMPORTS
#   define MEM_ALIGN(n) alignas(n)
#endif

#if (defined(WIN32) || defined(_WIN32) || defined(WINCE) || defined(__CYGWIN__))
#   define PACKED(__declare__) __pragma(pack(push,1)) __declare__ __pragma(pack(pop))
#else
#   define PACKED(__declare__) __declare__ __attribute__((__packed__))
#endif

#endif  //ARCSIM_API_HPP
