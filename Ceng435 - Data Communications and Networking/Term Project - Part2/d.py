ip_r1 = "10.10.3.2"            # ip address of the d-interface for receiving a message from r1
ip_r2 = "10.10.5.2"            # ip address of the d-interface for receiving a message from r2

ip_r1_send = "10.10.2.1"        # ip address of the interface for sending a msg to r1
ip_r2_send = "10.10.4.1"        # ip address of the interface for sending a msg to r2


port = 12346
chunk_size = 636                # To receive the messages in chunks of 636 bytes (header + payload)


import operator
import socket
import hashlib
from threading import Thread
from datetime import datetime
from threading import Lock
import time
import sys


last_sent_ack = 0        	# index of received_messages
packet_number = 0
received_messages = []			# Store the messages received from s to forward to d

chunks = [None] * 10000			# Store the ordered chunks

# Locks for avoiding race conditions
buffer_lock = Lock()


# Compute checksum value to compare with the one in the header
# and to compute checksum of the ACKs that will be sent.
# Returns 32 bytes hash
def checksum(message):
	m = hashlib.md5(message)
	return m.hexdigest()


# Packetize and create ACK message.
def createACK(packet_number):
	message = "ACK"
	
	packnum_str = str(packet_number)
	while len(packnum_str) < 4:			# Since chunks are at most 600 bytes, packet number can have at most 4 digits
										# If less, pad it with zeros
		packnum_str = "0" + packnum_str

	packet = packnum_str + message 		# Add packet number to the beginning of the message

	chk = str(checksum(packet))			# Compute checksum

	packet = chk + packet 				# Add checksum value to the beginning of the message

	return packet


# Checks if the checksum value in the header matches the checksum value
# computed by the d node.
# Returns true if the received packet is corrupted.
def isCorrupted(message, checksum_header):
	if str(checksum(message)) != str(checksum_header):
		return True
	return False


# Receives messages from r1 and r2
def receiveMsgFromRouters(sock):
	global packet_number, chunks

	while True:     
		
		received_message, address = sock.recvfrom(chunk_size)				# received_message stores the data received
																			# address stores the (ip, port) tuple of the sender(r1 or r2)
		if received_message:
			pack_num = int(received_message[32:36])
			#print "Packet", packet_num, "is received."

			buffer_lock.acquire()
			received_messages.append(received_message)						# Store received_message in received_messages
			chunks[pack_num] = received_message[36:]						# Keep data ordered
			sys.stdout.write(chunks[pack_num])
			buffer_lock.release()
		
			#if pack_num == 8333:
				#print time.time()

			packet_number += 1




# Sends the ACK message.
# Exploits both routers alternatingly (e.g. if previous ACK packet is sent to r1, sent this ACK packet to r2)
def sendMsgToRouters(sockdr1, sockdr2):
	global last_sent_ack

	while True:
		buffer_lock.acquire()
		while (last_sent_ack < len(received_messages)):			# While there is data to send
			received_message = received_messages[last_sent_ack]

			received_check = received_message[0:32]				# Checksum in the header
			pack_num = int(received_message[32:36])				# Packet number in the header

			if not isCorrupted(received_message[32:], received_check):
				if last_sent_ack % 2 == 1:           			# If last packet is sent to r1
					# Directly send this packet to the next hop(r2)
					sockdr2.sendto(createACK(pack_num), (ip_r2_send, port))
					#print "ACK SENT ! to r2"
				else:                                   		# If last packet is sent to r2
					# Directly send this packet to the next hop(r1)
					sockdr1.sendto(createACK(pack_num), (ip_r1_send, port))
					#print "ACK SENT ! to r1"
				last_sent_ack += 1
			#else:
				#print "C"

		buffer_lock.release()
	

# Creates sockets for sending and receiving messages.
# Starts the threads for sending and receiving messages.
def connection():
	sockr1d = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)             # Create a UDP socket for receiving messages from r1
	sockr1d.bind((ip_r1, port))                                            # Bind the socket to r1's address

	sockr2d = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)             # Create a UDP socket for receiving messages from r2
	sockr2d.bind((ip_r2, port))                                            # Bind the socket to r2's address

	sockdr1 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)             # Create a UDP socket for sending messages to r1

	sockdr2 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)             # Create a UDP socket for sending messages to r2

	receive_r1_thread = Thread(target = receiveMsgFromRouters, args = (sockr1d,))	# Create and start thread for receiving from r1
	receive_r1_thread.start()

	receive_r2_thread = Thread(target = receiveMsgFromRouters, args = (sockr2d,))	# Create and start thread for receiving from r2
	receive_r2_thread.start()

	send_thread = Thread(target = sendMsgToRouters, args = (sockdr1, sockdr2))		# Create and start thread for sending to r1 and r2
	send_thread.start()



if __name__ == "__main__":
	connection()
