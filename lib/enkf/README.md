### ENKF

ENKF(Integrated Model Inventory System) 是集合卡尔曼滤波库

### 编译

ENKF基于CMAKE进行编译，依赖CMAKE版本为3.x

```
mkdir build; cd build

FC=ifort cmake ../ # 普通编译
FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DTEST_FLAGS=test ../  # 单元测试

```

### 调试

```
gdb --args ./test/test.exe ../namelist.input
```

