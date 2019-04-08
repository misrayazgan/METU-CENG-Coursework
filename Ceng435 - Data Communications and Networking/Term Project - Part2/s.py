ip = "10.10.1.2"            # ip address of the next hop
ip_rcv = "10.10.1.1"        # ip address for receiving msg from b

port = 12346

chunk_size = 600           # To send the messages in chunks of 600 bytes

chunks = []

import socket
import struct
from threading import Thread
from datetime import datetime
import time


packet_number = 0


# Reads the file specified with its name by chunks.
# Chunk size is 600 bytes as specified above.
def readFile(filename):
	file = open(filename)
	
	while True:
		chunk = file.read(chunk_size)

		if chunk:
			chunks.append(chunk)
		else:
			break

			
# Target function for the sending thread.
# Sends the messages to the next hop address(broker).
def sendMsgToB(sock):
	global packet_number

	while packet_number < len(chunks):				# For each chunk in the file
		sock.send(chunks[packet_number])
		print "Packet", packet_number, "is sent."    # For debug
		packet_number += 1

		time.sleep(0.1)


# Creates the sockets for sending and receiving messages.
# Starts the threads for sending and receiving messages.
def connection():
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)    		# Create tcp socket
	sock.connect((ip, port))                                    		# Connect the socket to b's address (start the three-way handshake)


	send_thread = Thread(target = sendMsgToB, args = (sock,))			# Create thread for sending the message to b.
	send_thread.start()													# Start the sending thread.
	#print time.time() * 1000


if __name__ == "__main__":
	readFile("input.txt")						# Read the file in chunks
	connection()									# Start communication
