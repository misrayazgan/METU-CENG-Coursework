import itertools
from copy import deepcopy

# Predicate names, variables and function names start with a lower case letter.
# Constants starts with an upper case letter.

# While selecting clauses, set of support strategy is used:
# At least one of the parent clauses should be from the negated goal or its descendants.

class Predicate:
    def __init__(self, data):
        # data: string e.g. ~p(A,f(x))

        self.data = data
        if data[0] == '~':
            self.negation = True
            index = data.index('(')
            self.predicate_name = data[1:index]
            self.positive_data = data[1:]
        else:
            self.negation = False
            index = data.index('(')
            self.predicate_name = data[0:index]
            self.positive_data = data

        # number of constants, variables, functions
        self.elements = []
        # Parse predicate and get its elements
        self.parse()

    def parse(self):
        inside = self.positive_data[2:-1]
        stack = []
        split_indices = []
        for i, char in enumerate(inside):
            if char == '(':
                stack.append(char)
            elif char == ')':
                stack.pop()
            elif len(stack) == 0 and char == ',':
                split_indices.append(i + 2)

        # Get elements by splitting at split_indices
        start_index = 2
        for index in split_indices:
            element = self.positive_data[start_index: index]
            self.elements.append(element)
            start_index = index + 1

        self.elements.append(self.positive_data[start_index:-1])
        self.size = len(self.elements)

    def isEqual(self, other):
        if self.size != other.size:
            return False
        if self.data == other.data:
            return True

        return False

    def isNegated(self, other):
        if self.size != other.size:
            return False
        if self.positive_data == other.positive_data and self.negation != other.negation:
            return True

        return False

    # Return true if all elements of the predicate is constants
    def isConstantArgs(self):
        for element in self.elements:
            if isFunction(element):
                args = getFunctionArgs(element)
                for arg in args:
                    if not isConstant(arg):
                        return False
            elif not isConstant(element):
                return False

        return True

    # Return true if self is more general than other
    # Given predicates have the same size, name and negation
    # Compare their elements
    def isGeneralThan(self, other):
        for e1, e2 in zip(self.elements, other.elements):
            if isConstant(e1) and not isConstant(e2):
                return False

        return True

    # Replace occurences of old with new
    def replaceData(self, old, new):
        self.data = self.data.replace(old, new)
        self.positive_data = self.positive_data.replace(old, new)

    # Replace all data with respect to substitution
    def substituteData(self, substitution):
        if substitution != None:
            for element in self.elements:
                if element in substitution:
                    self.replaceData(element, substitution[element])
                    element = substitution[element]
        return self

    # Return disjunction of predicates in the predicate_list
    def disjunct(self, predicate_list):
        clause = ""
        for pred in predicate_list:
            clause += pred.data + ','

        clause = clause[:-1]
        return clause



class Clause:
    def __init__(self, data, parents = []):
        # data: string, one or more predicates separated by a comma
        # parents: list of Clause objects, two parent clauses for each resolvent clause

        self.data = data
        self.parents = parents
        self.predicates = []

        if self.data == "":
            # If empty clause, do not parse predicates
            self.data = "empty_clause"
        else:
            # Parse clause and initialize its predicates
            self.parse()

    def parse(self):
        stack = []
        split_indices = []

        # Find the indices of ','
        for i, char in enumerate(self.data):
            if len(stack) == 0 and char == ',':
                split_indices.append(i)

            if char == '(':
                stack.append(char)
            elif char == ')':
                stack.pop()

        # Get predicates by splitting at ','
        start_index = 0
        for index in split_indices:
            predicate = Predicate(self.data[start_index: index])
            self.predicates.append(predicate)
            start_index = index + 1

        self.predicates.append(Predicate(self.data[start_index:]))

    def isEqual(self, other):
        if self.data == other.data:
            return True
        return False

    # Return true if given predicate is in self.predicates
    def contains(self, pred):
        for p in self.predicates:
            if pred.isEqual(p):
                return True

        return False

    # Return true if self is a subset of the other
    def isSubset(self, other):
        for pred1 in self.predicates:
            if not other.contains(pred1):
                return False

        return True

    # Return true if clause is a tautology
    def isTautology(self):
        for pred1 in self.predicates:
            for pred2 in self.predicates:
                if pred1.isNegated(pred2):
                    return True

        return False

    # Print one line of resolution
    def getResolution(self):
        return self.parents[0].data + "$" + self.parents[1].data + "$" + self.data + "\n"

    # Return list of predicates with the given name and negation
    def getPredicatesByName(self, name, negation):
        return [pred for pred in self.predicates if pred.predicate_name == name and pred.negation == negation]

    # Apply substitution to all predicates of the clause
    def substituteData(self, substitution):
        if substitution != None:
            self.data = ""
            for pred in self.predicates:
                pred.substituteData(substitution)
                self.data += pred.data + ','

            self.data = self.data[:-1]
        return self

    def resolve(self, other):
        resolvents = []
        # Use copies for not substituting the data of the parent clauses
        self_copy = deepcopy(self)
        other_copy = deepcopy(other)

        for pred1 in self_copy.predicates:
            for pred2 in other_copy.predicates:
                if pred1.predicate_name == pred2.predicate_name and pred1.negation != pred2.negation:
                    substitution = unification(pred1.elements, pred2.elements, {})

                    # If a predicate in clause1 is the negated form of a predicate in clause2
                    # disjunct all predicates in two clauses except pred1 and pred2
                    if substitution != None:
                        preds_to_disjunct = combinePredicateLists(self_copy.predicates, other_copy.predicates, pred1, pred2)

                        for pred in preds_to_disjunct:
                            pred.substituteData(substitution)

                        disjuncted = []
                        for pred in preds_to_disjunct:
                            appendPredicate(pred, disjuncted)

                        clause_string = pred1.disjunct(disjuncted)
                        resolvent_clause = Clause(clause_string, parents=[self, other])
                        resolvents.append(resolvent_clause)

        return resolvents


knowledge_base_clauses = []
goal_clauses = []

def readInput(filename):
    with open(filename, "r") as input:
        # strip() gets rid of the newline character at the end
        line = input.readline().strip()
        numbers = line.split(" ")
        base_number = int(numbers[0])
        goal_number = int(numbers[1])

        # Read knowledge base clauses from input file and create Clause instances
        for i in range(base_number):
            clause = Clause(input.readline().strip())
            knowledge_base_clauses.append(clause)

        # Read goal clauses from input file and create Clause instances
        for i in range(goal_number):
            clause = Clause(input.readline().strip())
            goal_clauses.append(clause)


# Variables start with an uppercase letter
def isVariable(x):
    if isinstance(x, str) and x[0].islower() and '(' not in x:
        return True
    return False


# Constants start with a lowercase letter
def isConstant(x):
    if isinstance(x, str) and x[0].isupper():
        return True
    return False


# Functions start with a lowercase letter
def isFunction(x):
    if isinstance(x, str) and x[0].islower() and '(' in x:
        return True
    return False


# Return arguments of a given function
# Assuming that functions are not nested !!!!!!!!
def getFunctionArgs(f):
    args = []
    for i, char in enumerate(f):
        if char != '(' and char != ')' and char != ',' and i != 0:
            args.append(char)

    return args


def combinePredicateLists(list1, list2, self_pred, other_pred):
    result = list1 + list2
    to_be_removed = [self_pred, other_pred]

    for pred1, pred2 in itertools.combinations(result, 2):
        if pred1.isEqual(pred2):
            appendPredicate(pred1, to_be_removed)
        elif pred1.size == pred2.size and pred1.predicate_name == pred2.predicate_name and pred1.negation == pred2.negation:
            if not pred1.isConstantArgs() and not pred2.isConstantArgs():
                if pred1.isGeneralThan(pred2):
                    appendPredicate(pred2, to_be_removed)
                elif pred2.isGeneralThan(pred1):
                    appendPredicate(pred1, to_be_removed)

    for pred in to_be_removed:
        result.remove(pred)

    return result


# Return true if x occurs anywhere in y
def occurCheck(x, y):
    if x == y:
        return True
    elif isFunction(y):
        return x[0] == y[0] or occurCheck(x, y.getFunctionArgs())
    elif isinstance(y, list):
        for arg in y:
            if occurCheck(x, arg):
                return True
    return False


# substitution is a dictionary e.g. { 'a': 'x', 'b': 'y' }
# Returns a substitution
def unifyVariable(var, x, substitution):
    if var in substitution:
        val = substitution[var]
        return unification(val, x, substitution)
    elif x in substitution:
        val = substitution[x]
        return unification(var, val, substitution)
    elif occurCheck(var, x):
        return None
    else:
        substitution[var] = x
        return substitution


# x and y can be a variable, constant, list or a function
# At the beginning, x and y contain list of predicates of the parent clauses.
# Returns a substitution
def unification(x, y, substitution):
    if substitution == None:
        return None
    elif x == y:
        return substitution
    elif isVariable(x):
        return unifyVariable(x, y, substitution)
    elif isVariable(y):
        return unifyVariable(y, x, substitution)
    elif isFunction(x) and isFunction(y):
        return unification(getFunctionArgs(x), getFunctionArgs(y), unification(x[0], y[0], substitution))
    elif isinstance(x, list) and isinstance(y, list):
        return unification(x[1:], y[1:], unification(x[0], y[0], substitution))
    else:
        # If x and y are different constants
        return None


# If derivable, write "yes" to first line
# Else, write "no" to first line
def writeOutput(filename, status):
    with open(filename, "a") as output:
        output.write(status + "\n")


# Write all resolutions that contribute to the solution
def writeResolution(resolvent, filename):
    clauses = resolvent.parents
    solution = [resolvent.getResolution()]

    while len(clauses) > 0:
        last = clauses[-1]
        clauses = clauses[:-1]
        if last.parents != []:
            clauses += last.parents
            resolution = last.getResolution()
            solution.append(resolution)

    solution.reverse()

    with open(filename, "a") as output:
        for s in solution:
            output.write(s)


# Eliminate clauses containing tautologies
# e.g. p(a),~p(a),f(B,y)
def eliminateTautologiesinKB():
    to_be_removed = []

    for clause in knowledge_base_clauses:
        for pred1 in clause.predicates:
            for pred2 in clause.predicates:
                if pred1.isNegated(pred2):
                    appendClause(clause, to_be_removed)

    for clause in to_be_removed:
        knowledge_base_clauses.remove(clause)


# Eliminate all sentences that are more specific than an existing sentence in the KB.
# e.g. if p(x), eliminate p(A) and p(A),q(B)
def eliminateSubsumptionsinKB():
    to_be_removed = []

    for clause1 in knowledge_base_clauses:
        for clause2 in knowledge_base_clauses:
            if not clause1.isEqual(clause2):
                # Check if clause1 is a subset of clause2 under substitution
                if len(clause1.predicates) <= len(clause2.predicates):
                    substitution = {}
                    for pred1 in clause1.predicates:
                        preds = clause2.getPredicatesByName(pred1.predicate_name, pred1.negation)
                        for pred2 in preds:
                            if not pred1.isEqual(pred2):
                                substitution = unification(pred1.elements, pred2.elements, substitution)
                                if substitution != None:
                                    # Copies are used for not substituting the original clause
                                    copy_clause1 = deepcopy(clause1)
                                    copy_clause2 = deepcopy(clause2)
                                    copy_clause1.substituteData(substitution)

                                    if copy_clause1.isSubset(copy_clause2):
                                        appendClause(clause2, to_be_removed)


                elif len(clause1.predicates) > len(clause2.predicates):
                    substitution = {}
                    for pred1 in clause1.predicates:
                        preds = clause2.getPredicatesByName(pred1.predicate_name, pred1.negation)
                        for pred2 in preds:
                            if not pred1.isEqual(pred2):
                                substitution = unification(pred1.elements, pred2.elements, substitution)
                                if substitution != None:
                                    # Copies are used for not substituting the original clause
                                    copy_clause1 = deepcopy(clause1)
                                    copy_clause2 = deepcopy(clause2)
                                    copy_clause2.substituteData(substitution)

                                    if copy_clause2.isSubset(copy_clause1):
                                        appendClause(clause1, to_be_removed)

    for clause in to_be_removed:
        kb_and_goal_set.remove(clause)


# For all newly added clauses, remove subsumptions
def eliminateSubsumptions(kb_and_goal_set, goal_set, new_clause):
    to_be_removed = []

    for clause in kb_and_goal_set:
        # Check if clause1 is a subset of clause2 under substitution
        if len(clause.predicates) <= len(new_clause.predicates):
            substitution = {}
            for pred in clause.predicates:
                new_preds = new_clause.getPredicatesByName(pred.predicate_name, pred.negation)
                for new_pred in new_preds:
                    if not pred.isEqual(new_pred):
                        substitution = unification(pred.elements, new_pred.elements, substitution)
                        if substitution != None:
                            # Copies are used for not substituting the original clause
                            copy_clause1 = deepcopy(clause)
                            copy_clause2 = deepcopy(new_clause)
                            copy_clause1.substituteData(substitution)

                            if copy_clause1.isSubset(copy_clause2):
                                appendClause(new_clause, to_be_removed)


        elif len(clause.predicates) > len(new_clause.predicates):
            substitution = {}
            for pred in clause.predicates:
                new_preds = new_clause.getPredicatesByName(pred.predicate_name, pred.negation)
                for new_pred in new_preds:
                    if not pred.isEqual(new_pred):
                        substitution = unification(pred.elements, new_pred.elements, substitution)
                        if substitution != None:
                            # Copies are used for not substituting the original clause
                            copy_clause1 = deepcopy(clause)
                            copy_clause2 = deepcopy(new_clause)
                            copy_clause2.substituteData(substitution)

                            if copy_clause2.isSubset(copy_clause1):
                                appendClause(clause, to_be_removed)

    for clause in to_be_removed:
        kb_and_goal_set.remove(clause)
        goal_set.remove(clause)


# Append clause to list if it is not already in it
def appendClause(clause, list):
    for c in list:
        if c.isEqual(clause):
            return False

    list.append(clause)
    return True


# Append predicate to list if it is not already in it
def appendPredicate(pred, list):
    for p in list:
        if p.isEqual(pred):
            return False

    list.append(pred)
    return True


def resolution(filename):
    eliminateTautologiesinKB()
    eliminateSubsumptionsinKB()

    # set of support = goal clauses
    kb_and_goal_set = knowledge_base_clauses + goal_clauses
    goal_set = goal_clauses

    # For each clause in the sets, resolve their predicates
    for i, goal in enumerate(goal_set):
        for j, kb_goal in enumerate(kb_and_goal_set):
            resolvents = goal.resolve(kb_goal)

            if len(resolvents) == 0:
                continue

            for resolvent_clause in resolvents:
                if resolvent_clause.data == "empty_clause":
                    writeOutput(filename, "yes")
                    writeResolution(resolvent_clause, filename)
                    return

                if not resolvent_clause.isTautology():
                    appendClause(resolvent_clause, kb_and_goal_set)
                    appendClause(resolvent_clause, goal_set)

                    eliminateSubsumptions(kb_and_goal_set, goal_set, resolvent_clause)
                #print "resolution is:", resolvent_clause.getResolution()

    writeOutput(filename, "no")


if __name__ == "__main__":
    readInput("input.txt")
    output = open("output.txt","w+")
    resolution("output.txt")
