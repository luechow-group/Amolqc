## Coding guidelines for the program "amolqc"

Any new code should adhere to the guidelines below.

### 1. Indentation and Language

Fortran free-format style is used. For readability of 
the code, indentation should be used without exception.
Make sure that 3 "spaces" are used for indentation, not "tabs"
(i.e. the "tab" key indents with "spaces", this is usually an 
editor setting")

Language of variables and all comments is **English**

### 2. Naming conventions
* Modules: lower camel case with "_m" suffix ("randomWalker_m"). 
File name and module name should be identical.
* Variables: lower camel case ("randomWalker") unless uppercase is required to follow mathematical formulae.
* Procedures: upper camel case ("OptimizeParameters()")
* Type Names: upper camel case with "_t" suffix ("RandomWalker_t")

#### 2.1 In Modules that define "classes"

(a class defines a type and "methods" or functions work on this type):

Module name and type name should coincide (but module
with "_m" appended and type with "_t")

all functions have "this" of type "class type" as first argument.

all functions have the name: "type_methodx(this,...)"
(to avoid extensive renaming with use or complicated
construction of generic functions)


### 3. program structure

#### 3.1 Module concept

Use object-oriented approach where possible, meaning
group data and data accessing routines in one module

Allow all F90 parameter type checking, and F90 array
support by adding an interface module to all F77 or
module free F90 routines. Denote the interface module
as "nameinterface" and put it into a separate file
denoted here name_i.f\[90]. This module file produced from
name_i.f is equivalent to a compiled header file in C++.
The object file is basically empty and need not be 
included in the linking process.

This interface is only used if all calls to free routines
have a preceding "use ...interface".

The construction of interface files also helps avoiding
circular dependencies in Makefiles (although circular
dependencies hint to a less than optimal module design) 

How to write new code for amolqc:
Think about data structures for the problem, and routines
using this structures. Together they form a new module
(or more than one if a separation is possible: build small
modules!)

Those data that are used by more than one routine should 
become the data part of the new module, those used by only
one should be local variables of that routine (possibly
using 'save' to keep contents for next call)
Common data are at least parameters of the problem. They
are set with input routines, and accessed with output
routines.

IMPORTANT: use
```Fortran
private

public :: name_of_all_public_routines
```
to prevent access to the module data!

use access functions: getData(), setData(x)

possible exception: data access is part of time-determining
  code. There direct access is sensible, but should be used
  with care.

Naming convention for global data of a module:

Use "m"

use public data only if absolutely necessary (possibly
for speed reasons, but consider compiler inlining)

This design works only with a hierarchical module dependency
(i.e. a dependency tree of 'use' statements).
The general design of the code should allow this dependency.
If this is not possible (i.e. a circular dependency of 
modules in unavoidable), the data part of the module
might be factored out as global data and be accessed from
different modules. A pure data module  may only depend on
other data modules, and a hierarchical layout of data is
always possible.

Use the interface approach (free functions, subroutines and
an additional file with 'interface's) should be used only
for existing code.

Every Module **must** contain "implicit none".

Using other modules:
- use as little module as necessary
- write always an "only" list. This helps identifying where an undeclared variable/function comes from
  and restricts access to the absolute necessary items
  (this is not necessary when the "used" module has no public data, and a well designed public interface)

Structure of a module:
```Fortran
module moduleName
use kinds_m, only: r8
use module1, only: list_of_used_routines_or_data
use module2, only: list_of_used_routines_or_data

implicit none
private

public :: list_of_public_routines_or_data

    ! here is the data

contains

    ! here are the routines

end module moduleName
```

#### 3.2. Guidelines for individual routines

Good documentation of the routine interface (i.e. arguments, return values)
is imperative.

Use program language for documentation as much as possible,
i.e.
* all arguments declared with correct "intent" (and brief one line explanation after "!")
* all restrictions on the arguments and/or module variables enforced with
"call assert(condition,"failure note")" directly after declaration part

Use assert also in the routine to ensure that a condition
for the following part of the routine is fulfilled. E.g. if some
lines require variable "a" to be strictly positive use
```
  if (asserts) call assert(a > 0,'routine_name: condition a > 0 not satisfied')
```
"call assert" has no performance penalty due to the global variable asserts.

### 4. Error handling

Never use "print*" or "write(*,*)" because it is a parallel program!
write only to "iul" (the open file unit for logging)

Never call "stop" or "exit". It is a parallel program!
Delegate aborting the program to the function "abortp".
This one function can be adapted to handle ALL requests
to terminate the program.

### 5. Testing

#### 5.1. Module test

Provide testing for all functions of a module using pFUnit - [example](utils/generic/tests/blockAllocator_tm.pf).
To add all \*.pf files in a directory as a test, put
```cmake
add_pFUnit_test("ModuleName" "${libraries}")
```
into CMakeLists.txt

#### 5.2. Overall test of functionality

After the new functionality of the program works and is tested,
provide a test input (.in and .wf) with a few seconds runtime
that uses that new functionality and **[add it to the testsuite](testsuite/README.md)**
