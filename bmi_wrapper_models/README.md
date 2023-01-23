## C++ based BMI for C and C++ models
Example illustrating C++ based BMI for models written in C and C++. The BMI code mostly stays the same with a few exceptions, which are highlighted
inside the code. 
The example finds a root of cubic equations (coefficients are prescribed through a config file or can be set via BMI SetValue function) using
the Newton-Raphson algorithm.

#### Build:
  - mkdir build && cd build
  - cmake ../ -DMODELC:BOOL=ON (for c-based models)
  - cmake ../ -DMODELCXX:BOOL=ON (for cxx-based models)
  - make

#### Run:
  - run `./build/model_c configs/config_file.txt` (for c-based models)
  - run `./build/model_cxx configs/config_file.txt` (for cxx-based models)
  
### Additional notes
  - the BMI implementation for C and CXX models is not sigficantly different. 
  - Main differences are in the Initialize and Update method, and notes (numbered as note #1, note #2, ...) inside the code are provided to highlight the places where changes are needed/done

  
