# README #

Feynman diagram color factor for SU(N) gauge group numerical evaluation software.

### Installation: ###

Install Eigen library [http://eigen.tuxfamily.org]
And specify path to it when configure if it is not standard

```
$ ./configure CPPFLAGS=-I/usr/local/include/eigen3
$ make
```

### Usage: ###

Specify number of colors as first argument

`$./StupidPainter <NC>`

than using standard input stream specify color structure using:

* f(1,2,3) - for structure constants

* [T(1)T(2)T(3)] - for trace of generators

* To improve readability allowed spaces and separators as `*`, `:`,
`;`

To run tests use
`source test.sh`
