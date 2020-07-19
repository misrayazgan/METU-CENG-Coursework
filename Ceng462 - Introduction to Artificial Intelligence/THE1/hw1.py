import copy
import numpy as np

MAX_MOVE_DEPTH = 31
current_move_depth = 0
goal_state_positions = {}
solution_path = []
solution_path_ida = []


class State:
	def __init__(self):
		self.grid = np.array([[0,0,0],[0,0,0],[0,0,0]])
		self.g = 0
		self.h = 0
		self.parent = None
		self.blank_position = [0, 0]

	# Calculate number of nodes traversed from the initial(start) node.
	def set_g_score(self):
		if self.parent != None:
			self.g = self.parent.g + 1

	# Manhattan Distance Heuristic Function:
	# Calculate sum of Manhattan distance of the tiles except the blank
	# to their positions in the goal state.
	def set_h_score(self):
		global goal_state_positions
		total_distance = 0

		for row in range(3):
			for col in range(3):
				tile_value = self.grid[row][col]

				if tile_value != 0:
					goal_row = goal_state_positions[tile_value][0]
					goal_col = goal_state_positions[tile_value][1]

					row_distance = abs(row - goal_row)
					col_distance = abs(col - goal_col)
					total_distance += (row_distance + col_distance)

		self.h = total_distance

	def set_f_score(self):
		self.f = self.h + self.g

	def calculate_blank_position(self):
		for i in range(3):
			for j in range(3):
				if self.grid[i][j] == 0:
					return [i, j]


	def set_blank_position(self):
		self.blank_position = self.calculate_blank_position()

	def is_equal(self, state):
		if np.all(self.grid == state.grid):
			return True
		return False


	# Return the index in the list if there is a state object with self.grid in the list
	# If the object is not in the list, return -1
	def find_list_index(self, l):
		for i, element in enumerate(l):
			if np.all(self.grid == element.grid):
				return i
		return -1


# Fail if move count exceeds MAX_MOVE_DEPTH(31)
def is_fail():
	global current_move_depth

	if current_move_depth >= MAX_MOVE_DEPTH:
		return True
	return False


def find_inversion_count(flat_array):
	count = 0

	for i in range(8):
		for j in range(i + 1, 9):
			if flat_array[i] != 0 and flat_array[j] != 0 and flat_array[i] > flat_array[j]:
				count += 1

	return count


def is_even(a):
	if a % 2 == 0:
		return True
	return False


# Find whether the start and goal states are even or odd
# If even-even, odd-odd then solvable
def is_solvable(start_state, goal_state):
	start = start_state.grid.flatten()
	goal = goal_state.grid.flatten()

	inversion_start = find_inversion_count(start)
	inversion_goal = find_inversion_count(goal)

	if is_even(inversion_goal) and is_even(inversion_start):
		return True
	elif not is_even(inversion_goal) and not is_even(inversion_start):
		return True
	else:
		return False


# Write grid of a state to the specified output file
def write_output(filename, state):
	with open(filename, "a") as output:
		for i in range(3):
			string = ""
			for j in range(3):
				if j != 2:
					string += (str(state.grid[i][j]) + " ")
				else:
					string += (str(state.grid[i][j]) + "\n")
			output.write(string)
		output.write("\n")


# If the puzzle is not solvable, "fail" should be written to the specified output file
def write_fail(filename):
	with open(filename, "a") as output:
		output.write("fail")


def print_path():
	global solution_path
	solution_path.reverse()

	write_output("outputA.txt", start_state)

	for state in solution_path:
		write_output("outputA.txt", state)


# Write solution path for A* to the output file
def get_solution_path(state):
	global solution_path
	if state.parent == None:
		return print_path()
	else:
		solution_path.append(state)
		return get_solution_path(state.parent)


# Reads the problem from teh input file
# Creates and configures the start state and the goal state
def read_puzzle(filename, start_state, goal_state):
	global goal_state_positions

	with open(filename, "r") as input:
		file = input.read()
		lines = file.split("\n")

		# If last lines are empty
	 	while len(lines[-1]) == 0:
			lines = lines[:-1]

		# Read and save initial state of the puzzle
		for i in range(3):
			row_array = lines[i].split(' ')

			for j in range(3):
				start_state.grid[i][j] = int(row_array[j])

		# Read and save goal state of the puzzle
		for i in range(3):
			row_array = lines[i + 3 + 1].split(' ')

			for j in range(3):
				goal_state.grid[i][j] = int(row_array[j])
				goal_state_positions[int(row_array[j])] = [i, j]

		# Configure goal_state
		goal_state.set_g_score()
		goal_state.set_h_score()
		goal_state.set_f_score()
		goal_state.set_blank_position()

		# Configure start_state
		start_state.set_g_score()
		start_state.set_h_score()
		start_state.set_f_score()
		start_state.set_blank_position()


# Find the positions that are switchable with the blank tile(0).
# Consider right, left, up and down according to blank_position.
def find_positions_to_move(blank_pos):
	move_to = []
	blank_row = blank_pos[0]
	blank_col = blank_pos[1]

	# If not on the bottom line
	if blank_row < 2:
		move_to.append([blank_row + 1, blank_col])
	# If not on the top line
	if blank_row > 0:
		move_to.append([blank_row - 1, blank_col])
	# If not on the right line
	if blank_col < 2:
		move_to.append([blank_row, blank_col + 1])
	# If not on the left line
	if blank_col > 0:
		move_to.append([blank_row, blank_col - 1])

	return move_to


# Swap the blank tile with a specified position
# Returns the new state, sets its grid and blank position
def swap_tiles(state, prev_position, next_position):
	temp_state = State()

	temp_grid = copy.deepcopy(state.grid)
	temp = temp_grid[prev_position[0]][prev_position[1]]
	temp_grid[prev_position[0]][prev_position[1]] = temp_grid[next_position[0]][next_position[1]]
	temp_grid[next_position[0]][next_position[1]] = temp

	temp_state.grid = temp_grid
	temp_state.blank_position = next_position

	return temp_state


# Find successors of a state
def find_successors(state):
	blank_pos = state.blank_position
	move_to = find_positions_to_move(blank_pos)
	successors = []

	for pos in move_to:
		successor = swap_tiles(state, blank_pos, pos)
		successor.parent = state
		successor.set_g_score()
		successor.set_h_score()
		successor.set_f_score()

		successors.append(successor)

	return successors


# Find the state with the min f_score
def find_min_f_index(states):
	f_scores = []

	for state in states:
		f_scores.append(state.f)

	return f_scores.index(min(f_scores))


def a_star():
	global current_move_depth

	# After expanding the current state, put it into closed list
	closed_list = []

	# The states that can be generated from current state are added into open list
	open_list = [start_state]

	if not is_solvable(start_state, goal_state):
		write_fail("outputA.txt")
		return

	if len(open_list) == 0:
		write_fail("outputA.txt")
		return

	while len(open_list) > 0:
		min_f_index = find_min_f_index(open_list)
		min_f_state = open_list.pop(min_f_index)
		closed_list.append(min_f_state)

		current_move_depth = min_f_state.g

		successors = find_successors(min_f_state)

		if np.all(min_f_state.grid == goal_state.grid):
			get_solution_path(min_f_state)
			return

		if len(successors) == 0:
			write_fail("outputA.txt")
			return
		else:
			for successor in successors:
				if np.all(successor.grid == goal_state.grid):
					get_solution_path(successor)
					return

				open_index = successor.find_list_index(open_list)
				closed_index = successor.find_list_index(closed_list)

				if open_index == -1 and closed_index == -1:
					open_list.append(successor)

				elif open_index != -1:

					if successor.f < open_list[open_index].f:
						open_list[open_index].f = successor.f
						open_list[open_index].g = successor.g
						open_list[open_index].h = successor.h
						open_list[open_index].parent = successor.parent

				elif closed_index != -1:

					if successor.f <= closed_list[closed_index].f:
						closed_list.pop(closed_index)
						open_list.append(successor)

		if is_fail():
			write_fail("outputA.txt")
			return


def print_ida_solution():
	global solution_path_ida

	for state in solution_path_ida:
		write_output("outputIDA.txt", state)


def dfs(g, f_bound):
	global solution_path_ida

	current_state = solution_path_ida[-1]
	f = g + current_state.h

	if np.all(current_state.grid == goal_state.grid):
		return "OK"

	if f > f_bound:
		return f

	min_value = float('inf')

	successors = find_successors(current_state)

	for successor in successors:
		if successor not in solution_path_ida:
			solution_path_ida.append(successor)
			result = dfs(g + 1, f_bound)

			if result == "OK":
				return "OK"

			if result < min_value:
				min_value = result

			solution_path_ida = solution_path_ida[:-1]

	return min_value


def ida_star():
	global solution_path_ida

	if not is_solvable(start_state, goal_state):
		write_fail("outputIDA.txt")
		return

	f_bound = start_state.f

	solution_path_ida.append(start_state)

	while True:
		result = dfs(start_state.g, f_bound)

		if result == "OK":
			print_ida_solution()
			return

		if result == float('inf'):
			write_fail("outputIDA.txt")
			return

		f_bound = result



if __name__ == "__main__":
	# Start and goal states are configured when the puzzles are read from input file.
	start_state = State()
	goal_state = State()

	# Read the initial state and the goal state of the puzzle.
	read_puzzle("input.txt", start_state, goal_state)

	# Create the output files to write.
	outputA = open("outputA.txt","w+")
	outputIDA = open("outputIDA.txt", "w+")

	a_star()
	ida_star()
