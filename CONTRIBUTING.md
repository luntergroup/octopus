Our hope is that octopus will become a community driven project, so please feel tree to contribute! 

Before contributing please read the documentation, and ask questions if you're unsure of anything.

##Style

Please try to adhere to the existing code style, in particular:

####Naming

* Variable names are all lower case and seperated_with_underscores
* Enum (scoped and unscoped) members are lowerCamelCase
* constexpr variables are lowerCamelCase
* Global variables are lowerCamelCase
* Type names (including template types) are UpperCamelCase
* Function names are all lower case and seperated_with_underscores

####Brackets

* Classes and functions start with an open bracket on the next line, e.g:

    class Foo
    {
        // code
    };
    
    int bar(int x)
    {
        return 0;
    }

* Loops and if statements open the bracket on the same line.

####Namespaces

* All code goes in the octopus namespace

####Misc

* Use 4 spaces for tab indentation.
