# Carbon Gain vs Hydraulic Risk Stomatal Optimization Model V 2.0

__Coding language:__ C++

__Authors:__ German Vargas G. & William R.L. Anderegg

__Contact:__ german.vargas@utah.edu

------------

## Coding conventions:

The more people use this model, the better is our representation of physiological processes that determine how plant hydraulics affect gas exchange. For this reason we want our code to follow certain guidelines that we all can understand and use to make the source code more readable and easier to use and modify.

------------
__2)__ Avoid name spaces, always do std::cout.
```{Cpp}
#include <iostream>

namespace std; // WRONG
int main{
    std::cout << "Correct! Yaayyy" << std::endl;
}
```
__3)__ Minimize the use of global variables, if possible define them in the main program or locally within functions.
__4)__ Unless is necessary to return multiple values from a function, always use the Pass-by-const-reference when declaring arguments.
__5)__ Keep the model modular. If adding a function that quantifies a process related to carbon assimilation add it within the module 05CAssimilation.