## pyWATVina Installation：

**pyWATVina installation steps** 

1. **Unzip the file pywatvina.zip**

`unzip pywatvina.zip`

2. **decompress the static library of pywatvina**

`cd watvina`

`ar -x libwatvina.a`

some `.o` format files are released in this step

3. **Recompile the files to a share library file which could be imported by python，python**

**With your own environment under the specific version of python and boost **

`g++ -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fwrapv -O2 -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fwrapv -O2 -g -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 *.o -o _watvina_wrapper.so -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options`

A  `_watvina_wrapper.so` dynamic library file.

4. **remove the .a .o files** 

`$rm *.a *.o`

5. **copy the`watvina`folder to your python lib path**

like `${your_python_envs_lib/python3.x/dist-pakcages/}`