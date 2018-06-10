==========================================================================
 CTTG5.2 : Continuous time Tremblay Group
==========================================================================

:Authors: Charles-David Hébert, Maxime Charlebois, Patrick Sémon 
:Date: $Date: 2018-08-10 $
:Revision: $Revision: 5.2.0 $
:Description: Description

Naming conventions
-------------------
See this site:
    https://docs.microsoft.com/en-us/dotnet/standard/design-guidelines/naming-guidelines
 

General guidelines before pull requests
----------------------------------------

Build
^^^^^^^^^^^^^^^^^^^^^^
* Must build with clang for serial mode, no warnings, except for external libraries
* Must build with g++ for MPI mode, no warings, except for external libraries
* Must pass all tests
* Must repredouce the results in "test/Simulations"


Naming conventions
^^^^^^^^^^^^^^^^^^^
* Must respect the naming conventions


Formatting
^^^^^^^^^^^^^^^^
* Formatted according to Visual Studio standard, for exemple in visual studio code, put the folling preferences:
* "C_Cpp.clang_format_style": "Visual Studio"
* "editor.formatOnSave": true



    