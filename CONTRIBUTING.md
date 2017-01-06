# How to Contribute / Conventions

Contributions from the public are welcome.  If you plan on adding a new feature
I ask that you open a pull request with "WIP" (work in progress) in the title
to let me and the rest of the world know that you are working on this feature.
For all new features/classes etc. I ask that your PR also include a test case
that will ensure that I know when your feature breaks.

## Code Formating

Someday, when I'm not lazy, I'll get this more automated, until then *try* to
be consistent with camel case for classes/types and snake_case for
variables/functions/members.  Private/protected members should end in
underscores and code that is not intended to see the light of day (*i.e.*
implementation details) should be within the "detail" namespace.  Also please
try to keep lines equal to or less than 80 characters.  I code on a small
screen most of the time and it messes up the word wrap if you exceed roughly
that many characters...
