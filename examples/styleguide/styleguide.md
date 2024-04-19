# Style Guide

## Background 
Below outlines our style guide, the set of conventions and rules for our code. Some of the guidelines are coding best practices (and thus prevent bugs), others are more stylistic. Regardless, all guidelines must be followed. 

## General

### Collaboration
Remember that other people use this codebase too! If you are changing code that you are unfamiliar with, please consult with a knowledgeable member of the team first. Often, this can be done by looking at the most recent author of a commit.

### Paths
Never use an user-specific path in your code (e.g., `Users/author/downloads/data/mask.nii`), as this means your code will not work on someone else's machine without modification. Instead, make these paths an argument to your script or function.

### Dependencies
If your code requires external packages or dependencies, make sure to explicitly declare these in the language-appropriate dependency file (e.g., in Python, this is done using a `requirements.txt` or `Pipfile`.) Without this, it is difficult for someone else to run your code, as they need to download the correct packages *and* package versions.

### Abbreviations
Abbreviations and acronyms in variables can make code difficult to follow and hence are a source of bugs. Please be as descriptive as possible when creating variable names.

### TODOs
All TODOs need to have an author and a condition for when the TODO can be removed.

`TODO(ruben): Convert to version 11 of this library when it is stable`

### Error handling

#### Messsages 
All error messages should go to [STDERR](https://en.wikipedia.org/wiki/Standard_streams). This distinguishes normal diagnostic logging from critical issues.

#### Exiting
If a critical issue occurs such that program execution must be halted, exit with a non-zero code and log a helpful error message. For example, in Bash, this can be done as the following:

```bash
if [[test -f "file.csv"]]; then
  echo "Fatal: file.csv does not exist" >&2
  exit 1
fi
```

### Comments
Comments should focus on the "why" rather than the "what." That is, the comment explains the motivation or necessity for code, rather than what the code is doing. As a general principle, code should be written so simply that the average programmer does not have trouble understanding it.

The exception to this rule is when comments are used as "bookmarks" in long block of code and help the reader navigate logic.

#### Code
Never commit commented-out code to `main`. Commented-out code makes the codebase hard to maintain, as it is not clear why the code was commented out (as opposed to removed) and whether the code is still usable. If you are concerned about losing code, remember that you can commit the code and push it to your remote feature branch. The code will still be available, as long as the branch is not deleted.

### Global variables
Mutable global variables are convenient, but often bug-prone. Use sparingly.

### Language/framework choice
The right tool for a problem will depend on the problem as well as the author's familiarity with the tool. However, when possible, prefer to use an open source language (e.g., Python) over a commercial one (e.g., Matlab) as open-source projects tend to have a better ecosystem of extensions. Also, open source is free!

### Printing
Be mindful of how much and what you are logging to STDOUT. While it is helpful to print events, variables, and outputs when developing, this should be done sparingly for finished code.
If it is necessary to print copious amounts of text, please consider providing a "verbose" option to your code, where `verbose=false` means that logs are reduced.

## Python

### Linting
Run `pylint` over your code using `pylint src`. `pylint` is a tool for finding bugs and style problems in Python source code. Some warnings may be spurious and these can be ignored.

### Indentation
Use tabs, not spaces.

## Bash/Shell scripts

### When to use Shell
Shell should only be used for small utilities or simple wrapper scripts. Any script over 100 lines should be broken up into smaller scripts, or moved to a structured language (e.g., Python).