# Heredity

This Python project is designed to calculate the probability of a person having a particular trait based on genetic inheritance and mutation rates. It takes a CSV file containing information about individuals, their parents, and their trait status, and computes the probability distributions for each individual having the trait and possessing the gene responsible for the trait.

## Usage

To run the program, execute the following command:
```txt
python heredity.py data.csv
```
Replace `data.csv` with the path to your input CSV file.

## Input File Format

The input CSV file should have the following columns:

- `name`: The name of the individual.
- `mother`: The name of the individual's mother (if known), or blank if unknown.
- `father`: The name of the individual's father (if known), or blank if unknown.
- `trait`: The trait status of the individual, either `1` (true) or `0` (false), or blank if unknown.

## Probability Calculations

The program uses a set of predefined probabilities for gene inheritance and trait expression:

- `PROBS["gene"]`: The unconditional probabilities of having 0, 1, or 2 copies of the gene.
- `PROBS["trait"]`: The probabilities of expressing the trait given 0, 1, or 2 copies of the gene.
- `PROBS["mutation"]`: The probability of a mutation occurring during gene inheritance.

The program calculates the joint probability of all possible combinations of individuals having the gene and the trait. It then updates the probability distributions for each individual based on these joint probabilities.

## Output

The program outputs the probability distributions for each individual, showing the probabilities of having the trait and possessing the gene (0, 1, or 2 copies).

## Dependencies

This project requires Python 3 and the `csv` and `itertools` modules from the Python standard library.
