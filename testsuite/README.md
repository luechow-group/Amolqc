### Adding Tests

For a new test you need a new folder in "run" directory contains a script named "RUN".
"RUN" can be written in any scripting language but make sure that you have the runtime environment
everywhere. For your script you may consider these points.

- It should be possible to generate reference file(if it is necessary) using "GEN" environment variable.

- Reference files should be saved with .$REF extension this helps to keep ref file for both parallel and non parallel
   calculations.

- You can use input template file also wave function to avoid unnecessary repeated files but it is not mandatory.

- in case of failure "FAILED=1" environment variable should be exported and a warning should be printed (red color code font)
   is preferred.

You can look at "run/01_sample" folder as an example.

### Updating References

**Only do this for a good reason! Has to be discussed with the Amolqc administrator!**

```
cd testsuite/run
chmod +x tests.sh
```
```
./tests.sh GEN
```
generates reference files for single core calculations

```
./tests.sh PARALLELGEN
```
generates reference files for calculations with two cores