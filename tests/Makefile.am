# Makefile.am


check_PROGRAMS = iedera_gcov

TESTS      = iedera_gcov_test1 iedera_gcov_test2 iedera_gcov_test3 iedera_gcov_test4 iedera_gcov_test5 iedera_gcov_test6
TESTS_LOGS = iedera_gcov_test1 iedera_gcov_test2 iedera_gcov_test3 iedera_gcov_test4 iedera_gcov_test5 iedera_gcov_test6

iedera_gcov_CXXFLAGS = -fprofile-arcs -ftest-coverage -g -O0
iedera_gcov_SOURCES = ../src/main.cc ../src/seed.cc ../src/macro.cc

EXTRA_DIST = iedera_gcov_test1 iedera_gcov_test2 iedera_gcov_test3 iedera_gcov_test4 iedera_gcov_test5 iedera_gcov_test6
