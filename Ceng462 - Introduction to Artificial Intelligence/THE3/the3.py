from copy import deepcopy
import operator

class Node:
    def __init__(self, type, id):
        self.type = type
        self.id = id


class Action:
    def __init__(self, reward, probs_dict):
        self.reward = reward
        # Dictionary in form {node_id: probability}
        self.probs_dict = probs_dict


class Transition:
    def __init__(self, source_id, dest_id, reward):
        self.source_id = source_id
        self.dest_id = dest_id
        self.reward = reward


round_nodes = []        # this universe
vortex_nodes = []       # cross nodes/blackholes
room_nodes = []         # square nodes/transition nodes
star_nodes = []         # other universe
goal_node = None        # Skynet

learning_rate = 0
discount_factor = 0
transitions = []
transition_rewards = {}        # Store {(source_id, dest_id): reward}
actions_from_nodes = {}     # Store {node_id: action list}
actions = {}                # Store {(action_id, node_id): (reward, prob list)}


def read_input(filename):
    global learning_rate, discount_factor, actions, transitions, transition_rewards, actions_from_nodes, goal_node

    with open(filename, "r") as input:
        # Read node types
        node_types = input.readline().strip()
        for id, type in enumerate(node_types):
            node = Node(type, id)
            if type == 'R':
                round_nodes.append(node)
            elif type == 'S':
                star_nodes.append(node)
            elif type == 'V':
                vortex_nodes.append(node)
            elif type == 'O':
                room_nodes.append(node)
            elif type == 'G':
                goal_node = node

        # Read learning rate and discount factor
        line = input.readline().strip().split(" ")
        learning_rate = float(line[0])
        discount_factor = float(line[1])

        # Read deterministic transitions
        number_of_transitions = int(input.readline().strip())

        for i in range(number_of_transitions):
            line = input.readline().strip()
            args = line.split(" ")
            transition = Transition(int(args[0]), int(args[1]), int(args[2]))
            transitions.append(transition)
            transition_rewards[(transition.source_id, transition.dest_id)] = transition.reward

        # Read actions for each node in other universe except Skynet(goal node)
        number_of_actions = int(input.readline().strip())

        for i in range(number_of_actions):
            line = input.readline().strip()
            args = line.split(" ")
            actions_from_nodes[int(args[0])] = [int(arg) for arg in args[1:]]

        # Read the effects of these actions
        line = input.readline().strip()

        action = None
        action_id = 0

        while line != "E":
            if line[0] == "a":
                # New action
                args = line.split(" : ")
                action_id = int(args[1])

                node_id = int(input.readline().strip())
                reward = int(input.readline().strip())

                probs_dict = {}
                line = input.readline().strip()
                while line != "$":
                    node, prob = line.split(" ")
                    line = input.readline().strip()
                    probs_dict[int(node)] = int(prob) / 100.0

                action = Action(reward, probs_dict)
                actions[(action_id, node_id)] = action

            elif line == "$":
                # New node in action
                line = input.readline().strip()
                if line == "#":
                    continue

                node_id = int(line)
                reward = int(input.readline().strip())

                probs_dict = {}
                line = input.readline().strip()
                while line != "$":
                    node, prob = line.split(" ")
                    line = input.readline().strip()
                    probs_dict[int(node)] = int(prob) / 100.0

                action = Action(reward, probs_dict)
                actions[(action_id, node_id)] = action

            elif line == "#":
                # Action ended
                line = input.readline().strip()


q_table = {}        # Store {(source_id, dest_id): q_value}
copy_q_table = {}

# Initialize Q table with 0s
def init_Qtable():
    global q_table, transition_rewards

    zeros = [0] * len(transition_rewards.keys())
    q_table = dict(zip(transition_rewards.keys(), zeros))
    copy_q_table.update(deepcopy(q_table))


# Find the max Q value for the next state
def find_next_max_q(next_id):
    global q_table, copy_q_table
    q_values = []

    # Find the nodes that can be reached from next node and store all possible Q values
    for source_id, dest_id in q_table.keys():
        if source_id == next_id:
            q_values.append(copy_q_table[(source_id, dest_id)])

    # Find max Q value
    if len(q_values) > 0:
        return max(0, max(q_values))
    else:
        return 0


# Episodes are the node_ids
def calculate_Qtable(node_ids):
    global q_table, transition_rewards, copy_q_table

    room_ids = [node.id for node in room_nodes]
    vortex_ids = [node.id for node in vortex_nodes]

    for i, current_id in enumerate(node_ids):
        # Room nodes are the goals for this part, so do not consider the last node id.
        # Vortex nodes are the blackholes.
        if i != len(node_ids) - 1 and not current_id in room_ids + vortex_ids:

            next_id = node_ids[i + 1]
            next_max = find_next_max_q(next_id)
            current_reward = transition_rewards[(current_id, next_id)]

            new_q = (1 - learning_rate) * q_table[(current_id, next_id)] + learning_rate * (current_reward + discount_factor * next_max)
            q_table[(current_id, next_id)] = new_q
        else:
            break

    copy_q_table = deepcopy(q_table)


def print_Qtable():
    global q_table, transition_rewards

    shape = (len(round_nodes) + len(room_nodes) + len(vortex_nodes), len(round_nodes) + len(room_nodes) + len(vortex_nodes))

    printed_q = []

    for i in range(shape[0]):
        line = ""
        for j in range(shape[1]):
            if not (i,j) in transition_rewards.keys():
                line += '_ '
            else:
                line += str(q_table[(i,j)]) + ' '

        print line


v_table = {}           # Store {node_id: value}

def init_Vtable():
    global v_table, room_nodes, star_nodes, goal_node

    all_nodes = room_nodes + star_nodes + [goal_node]
    all_node_ids = [node.id for node in all_nodes]

    zeros = [0] * len(all_nodes)
    v_table = dict(zip(all_node_ids, zeros))


# Implemented for other universe (stochastic case)
# Returns V table
def value_iteration(epsilon=0.001):
    global v_table

    nodes = room_nodes + star_nodes
    node_ids = [node.id for node in nodes]

    while True:
        delta = 0
        new_v = v_table.copy()

        # Excluding goal node
        for node_id in node_ids:
            values = []
            for action_id in actions_from_nodes[node_id]:
                action = actions[(action_id, node_id)]

                probvalue_sum = 0
                for node in action.probs_dict.keys():
                    prob = action.probs_dict[node]
                    probvalue_sum += prob * new_v[node]

                value = action.reward + discount_factor * probvalue_sum
                values.append(value)

            v_table[node_id] = max(values)

            if abs(v_table[node_id] - new_v[node_id]) > delta:
                delta = abs(v_table[node_id] - new_v[node_id])

        # Convergence check
        if delta < epsilon * (1 - discount_factor) / discount_factor:
            return


def find_expected_value(action_id, node_id):
    global v_table, actions

    action = actions[(action_id, node_id)]
    sum = 0

    for node in action.probs_dict.keys():
        prob = action.probs_dict[node]
        sum += prob * v_table[node]

    return sum


def find_policy():
    global v_table

    nodes = room_nodes + star_nodes
    node_ids = [node.id for node in nodes]

    # Store {node_id: best_action}
    policy = {}

    for node_id in node_ids:
        values = {}
        for action_id in actions_from_nodes[node_id]:
            value = find_expected_value(action_id, node_id)
            values[action_id] = value

        policy[node_id] = max(values.iteritems(), key=operator.itemgetter(1))[0]

    return policy


# policy contains the best nodes that can be reached from given node
def print_policy(policy):
    nodes = room_nodes + star_nodes
    node_ids = [node.id for node in nodes]

    for node_id in node_ids:
        to_be_printed = str(node_id) + " "
        action = actions[(policy[node_id], node_id)]

        for node in action.probs_dict:
            to_be_printed += str(node) + ','

        print to_be_printed[:-1]

    print "\n"


def print_Vtable():
    global v_table

    for node_id in v_table.keys():
        print node_id, v_table[node_id]

    print "\n"


def interactive_session():
    dollar_count = 0

    while True:
        input = raw_input()

        if input == '$':
            dollar_count += 1

        if dollar_count == 0:
            node_ids = input.split(" ")
            node_ids = [int(id) for id in node_ids]

            # Implemented for this universe (deterministic case)
            # Goal nodes are the rooms, do not consider further.
            calculate_Qtable(node_ids)
            print_Qtable()

        elif input == '$' and dollar_count == 1:
            value_iteration()
            policy = find_policy()
            print_Vtable()
            print_policy(policy)

        elif dollar_count == 1 and input == 'c':
            # Q learning is terminated
            # print value table and policy for first iteration
            value_iteration()
            policy = find_policy()
            print_Vtable()
            print_policy(policy)

        elif dollar_count == 2:
            break


if __name__ == "__main__":
    read_input("the3.inp")

    init_Qtable()
    init_Vtable()
    interactive_session()
