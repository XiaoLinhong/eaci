### EACI
Emission Adjustment Coefficient Inversion
### 编译

EACI 基于CMAKE进行编译，依赖CMAKE版本为3.x

```
mkdir build; cd build

FC=ifort cmake ../ # 普通编译
FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug ../ # 调试

make
```

### 调试

```
gdb --args ./eaci.exe ../namelist.input
```

### 执行

```
./eaci.exe ../namelist.input
```
