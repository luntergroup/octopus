Our hope is that octopus will become a community driven project, so please feel tree to contribute! 

## General guidance

Before contributing please read the documentation, check existing issues and branches. Please ask questions if you're unsure of anything.

We follow [this](http://nvie.com/posts/a-successful-git-branching-model/) git branch model. In particular:

* The `master` branch is for stable versions only. Usually only tagged releases will be commited to `master`.
* The `develop` branch is the main working branch. Only small changes should be commited directly to this branch.
* New major features should branch from `develop` and be named like `feature/name`.
* Minor bug fixes should branch from the relevant branch and be named like `fix/name`.
* Experimental work should branch from the relevant branch and be named like `exp/name`.

## Style

Please try to adhere to the existing code style, in particular:

#### Naming

* Variable names are all lower case and seperated with underscores (e.g. `foo_bar`).
* Type names (including template types and aliases) are UpperCamelCase (e.g. `FooBar`).
* Function names are all lower case and seperated_with_underscores (e.g. `foor_bar()`).
* `enum` and `enum class` members are all lower case and seperated with underscores (e.g. `foo_bar`).
* `constexpr` variables are lowerCamelCase (e.g. `fooBar`).
* Global variables are lowerCamelCase (e.g. `fooBar`).

#### Brackets

* Classes and functions start with an open bracket on the next line, e.g:

```cpp
class Foo
{
    // code
};

int bar(int x)
{
    return 0;
}
```

* Loops and if statements open the bracket on the same line, e.g:

```cpp
for (;;) {
    // stuff
}

if (something) {
    // stuff
} else {
    // other stuff
}
```

* Small `if` statements can appear on the same line if it improves readability:

```cpp
if (something) return true;
```

#### Namespaces

* All code goes in the `octopus namespace`.
* The open bracket appears on the same line as the `namespace` name (e.g. `namespace octopus {`).

#### Misc

* Use 4 spaces for tab indentation.
* There is no strict limit for line lengths; use best judgement. 
