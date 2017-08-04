AC_PREREQ([2.61])
AC_INIT([StupidPainter], [1.0], [pikelner@theor.jinr.ru])
AM_INIT_AUTOMAKE([foreign])

AC_PROG_MAKE_SET
AC_PROG_CXX
AC_LANG(C++)

CXXFLAGS="$CPPFLAGS -std=c++11 -fPIC -O3"

AC_CONFIG_MACRO_DIR([m4])

AX_CXX_COMPILE_STDCXX(11, noext, mandatory)

dnl Eigen
AC_CHECK_HEADERS([Eigen/Dense],
                 [],          
                 [AC_MSG_ERROR([Could not find the Eigen header files.])]
                )

dnl GNU Readline
AC_CHECK_HEADERS([readline/readline.h readline/history.h],
                 [],          
                 [AC_MSG_ERROR([Could not find the GNU readline header files.])]
                )
AC_CHECK_LIB([readline], [readline])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
