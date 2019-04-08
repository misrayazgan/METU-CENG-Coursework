import vulnerable
import sys

def pad_DP(string):
	return string + '\x01'


def findCipher(padded_DP):
	T = [0] * 16

	ciphertext = list(original_ciphertext)
	cipher_to_return = list(original_ciphertext)

	padding_length = 1;

	cipher_index = 31


	for j in range(0, len(padded_DP)):					# for each character of the 2. ciphertext block C[16:31]

		success_array = []
		cipher_index = 31 - j
		
		for i in range(0, 256):							# try all possible hex values 0-255

			ciphertext[cipher_index] = chr(i)			

			if vulnerable.decr(''.join(ciphertext)) == "SUCCESS":
				success_array.append(chr(i))
				if len(success_array) == 2:
					break
		

		if len(success_array) == 2:
			for m in range(2):
				if success_array[m] != original_ciphertext[cipher_index]:
					T[15 - j] = ord(success_array[m]) ^ padding_length				# T[15] = i xor 0x1, T[14] = i xor 0x2, ...
					padding_length += 1

					for k in range(0, j + 1):										# update ciphertext values: C[31] = T[15] xor 0x2, ...
						ciphertext[31 - k] = chr(T[15 - k] ^ padding_length)

		elif len(success_array) == 1:
			T[15 - j] = ord(success_array[0]) ^ padding_length						# T[15] = i xor 0x1, T[14] = i xor 0x2, ...
			padding_length += 1

			for k in range(0, j + 1):												# update ciphertext values: C[31] = T[15] xor 0x2, ...
				ciphertext[31 - k] = chr(T[15 - k] ^ padding_length)
		

	for i in range(0, len(padded_DP)):
		cipher_to_return[31 - i] = chr(T[15 - i] ^ ord(padded_DP[len(padded_DP) - 1 - i]))


	return ''.join(cipher_to_return)



if __name__ == "__main__":

	original_ciphertext = sys.argv[1]
	
	DP = sys.argv[2]
	
	padded_DP = pad_DP(DP)
	
	cipher = findCipher(padded_DP)

	file = open("cipher.txt", "w+")
	file.write(cipher)


