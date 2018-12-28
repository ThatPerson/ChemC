# ChemC

ChemC is an attempt at porting [ChemKit](https://www.github.com/ThatPerson/ChemKit) from Python to C. It can do most of the same things and is a little bit easier to use, but it shares many of the fatal flaws with ChemKit.

ChemC works as a command line program. It can take the following commands;

| Command | Example Arguments | Explanation                           |
|-|-|-|
| help| | Lists help. |
| element| C | Prints out a bunch of information, some from the periodic table and some predicted |
| length | O H | Outputs the sum of the radii of the two elements |
| mass | O H | Outputs the reduced mass |
| compound | H2O | Outputs a load of junk. In this junk there is the molecular mass, the programs assumptions at the bonding in the molecule (which is fairly okay for small molecules like H2O, is utterly awful for large molecules like C6H12O6), predicted force constants for the bonds it thinks exists, and then predicted melting points (using the same awful models as in ChemKit). Not really useful but it's sometimes interesting. |
| reaction | AlCl3 Fe | Attempts to predict the products of the given reaction. This isn't done using any mechanistic understanding. |

