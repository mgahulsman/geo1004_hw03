# GEO1004 - Homework 3 - BIM to Geo conversion using voxels

## Team members
Anna Scherrenburg - 5212723

Maarten Hulsman - 5379989

## Code instructions
Before running the script, ensure your input data is placed in the correct relative directory. The program expects an input file to be present in the ../input/ folder relative to the execution path. In order to write to the CityJSON file, it is important to have json.hpp in the include folder.

```text
.
├── input/
│   └── input_file.obj
├── src/
│   └── main.cpp
├── include/
│   └── json.hpp
├── CMakeLists.txt
└── README.md
```

As long as the input file is in the right place, you can execute the program directly without passing any additional arguments. The core logic relies on a central parameters section within the source code.

Grid Size: This is the most important variable in the pipeline. It controls the spatial resolution of the processing grid. To adjust the grid size, you must modify its value manually inside the code. Look for the configuration block at the top of your main script:
