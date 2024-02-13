import csv
import itertools
import sys

#probabilities involing the given genes
PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    #PS: variables representing probabilities will be initialized to 1 to ensure multiplication netrulity during first manipulations
    joint_p = 1

    #loop over people dict
    for person in people:
        #store information about perosn`s genetics`
        genes_prob = 1
        genes_in_person = (2 if person in two_genes else 1 if person in one_gene else 0)
        has_trait  = person in have_trait

        #info about parents of person
        person_father = people[person]["father"]
        person_mother = people[person]["mother"]

        #case: no listed parents => use standard probabilities
        if not person_mother and not person_father:
            genes_prob = PROBS['gene'][genes_in_person]
        #case: parents listed => calculte probabilities of genes from inheritance and parents information]
        else:
            #prob of each parent passing the considered gene
            passed_gene_father = inherit(one_gene, two_genes, person_father)
            passed_gene_mother = inherit(one_gene, two_genes, person_mother)

            #compute the probability for the gene
            if genes_in_person == 2:
                genes_prob *= passed_gene_father * passed_gene_mother
            elif genes_in_person == 1:
                genes_prob *= (1-passed_gene_father)*passed_gene_mother + (1-passed_gene_mother)*passed_gene_father
            #case: person without copies of the mutated gene
            else:
                genes_prob *= (1-passed_gene_father) * (1-passed_gene_mother) 

        #consider probabilities of person having and not having the trait, given the probabilities of the number of genes
        genes_prob *= PROBS['trait'][genes_in_person][has_trait]
        joint_p *= genes_prob

    return joint_p

def inherit(one_gene, two_gene, parent):
    """
    Based on the number of genes a parent has and the probability of mutation, returns the probability that a child will receive a copy of the harmful gene. 
    """
    passing_gene_prob = 0.0
    
    #case: parent has 2 mutated genes
    if parent in two_gene:
        passing_gene_prob = (1-PROBS['mutation'])
    #case: parent has 1 mutated gene 
    elif parent in one_gene:
        passing_gene_prob = 0.5
    #case: parent has no copy of mutated gene
    else:
        passing_gene_prob = PROBS['mutation']

    return passing_gene_prob

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    #iterate over people in people dict
    for person in probabilities:
        #collect informaiton about current person
        person_num_genes = (2 if person in two_genes else 1 if person in one_gene else 0)
        has_trait = person in have_trait

        #update probabilities for current person
        #update gene probability
        probabilities[person]['gene'][person_num_genes] += p
        probabilities[person]['trait'][has_trait] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    #iterate over probabilities dict
    for person in probabilities:
        #get total sum of trait and genes probabilities
        total_trait = sum(probabilities[person]['trait'].values())
        total_genes = sum(probabilities[person]['gene'].values())

        #normalize values by dividing by the sums above
        probabilities[person]['trait'] = {trait: (p/total_trait) for trait, p in probabilities[person]['trait'].items()}
        probabilities[person]['gene'] = {gene: (p/total_genes) for gene, p in probabilities[person]['gene'].items()}


if __name__ == "__main__":
    main()
