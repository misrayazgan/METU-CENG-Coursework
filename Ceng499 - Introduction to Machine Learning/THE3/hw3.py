import numpy as np 
import sys



def readFile(filename):
	with open(filename, "r") as inputfile:
		f = inputfile.read()
		lines = f.split("\n")
		
		lines = lines[:-1]			# last line is empty

		columns = lines[0].split()
		shape = (len(lines), len(columns))
		matrix = np.zeros(shape, float)

		for i, line in enumerate(lines):
			row = line.split(' ')
			for j, number in enumerate(row):
				matrix[i][j] = float(number)

		return matrix



def forward(observations, estimate, transition):
	T = observations.shape[1]		# number of observations
	N = transition.shape[0]			# number of states including start and end states

	forward_matrix = np.zeros((N, T), float)

	# For the first observation, calculate 
	# (prob. of moving from start state to other states) * (prob. of obs0 happening in that state)
	for s in range(1, N - 1):
		forward_matrix[s][0] = transition[0][s] * estimate[s][int(observations[0][0])]


	# For other observations(time steps)
	for t in range(1, T):
		for s in range(1, N - 1):
			for i in range(1, N - 1):
				forward_matrix[s][t] += forward_matrix[i][t - 1] * transition[i][s] * estimate[s][int(observations[0][t])]

	
	# For the last observation, calculate sum of
	# (prob. of moving from other states to end state) * (prob. of having the previous observation in other states)
	for s in range(1, N - 1):
		forward_matrix[N - 1][T - 1] += forward_matrix[s][T - 1] * transition[s][N - 1]

	return forward_matrix[N - 1][T - 1]





def viterbi(observations, estimate, transition):
	T = observations.shape[1]		# number of observations
	N = transition.shape[0]			# number of states including start and end states

	viterbi_matrix = np.zeros((N, T), float)
	back_pointer = np.zeros((N, T), int)

	# For the first observation, calculate
	# (prob. of moving from start state to other states) * (prob. of obs0 happening in that state)
	for s in range(1, N - 1):
		viterbi_matrix[s][0] = transition[0][s] * estimate[s][int(observations[0][0])]
		back_pointer[s][0] = -1			# Previous state is start state


	# For other observations(time steps)
	for t in range(1, T):
		for s in range(1, N - 1):		# states excluding start and end
			probs = []
			state_probs = []
			for i in range(1, N - 1):
				prob = viterbi_matrix[i][t - 1] * transition[i][s] * estimate[s][int(observations[0][t])]
				probs.append(prob)

				state_prob = viterbi_matrix[i][t - 1] * transition[i][s]
				state_probs.append(state_prob)

			viterbi_matrix[s][t] = max(probs)
			back_pointer[s][t] = state_probs.index(max(state_probs)) + 1


	# For the last observation
	probs = []
	state_probs = []
	for s in range(1, N - 1):
		prob = viterbi_matrix[s][T - 1] * transition[s][N - 1]
		probs.append(prob)

		state_prob = viterbi_matrix[s][T - 1] * transition[s][N - 1]
		state_probs.append(state_prob)

	viterbi_matrix[N - 1][T - 1] = max(probs)
	back_pointer[N - 1][T - 1] = state_probs.index(max(state_probs)) + 1


	# Get the best path
	best_states = []
	best_states.append(back_pointer[N - 1][T - 1])
	for t in range(T -1, 0, -1):
		previous_state = back_pointer[best_states[-1]][t]
		best_states.append(previous_state)

	return best_states[::-1]




if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "Exactly 4 arguments should be given!"
		sys.exit()

	transition = readFile(sys.argv[2])
	estimate = readFile(sys.argv[3])
	observations = readFile(sys.argv[4])

	if sys.argv[1] == "forward":
		# call forward algorithm
		print forward(observations, estimate, transition)
	elif sys.argv[1] == "viterbi":
		# call viterbi algorithm
		print viterbi(observations, estimate, transition)



