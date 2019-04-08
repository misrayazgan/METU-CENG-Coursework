ip_s = "10.10.1.2"          # ip address of the b-interface for receiving a message from s
ip_r1 = "10.10.3.2"         # ip address of the next hop(r1)
ip_r2 = "10.10.5.2"         # ip address of the next hop(r2)

ip_r1_rcv = "10.10.2.1"     # ip address of the interface for receiving a msg from r1
ip_r2_rcv = "10.10.4.1"     # ip address of the interface for receiving a msg from r2


port = 12346

chunk_size = 600             # To receive the messages in chunks of 600 bytes


import socket
import time
import hashlib
from threading import Thread
from threading import Lock
import sys
from datetime import datetime


packet_number = 0

send_buffer = []			# Store the messages received from s to forward to d
last_sent_packet = 0        # index of send_buffer
base = 0					# index of first nonACKed packet - 1
window_size = 10

received_acks = [0] * 10000			# Store the ACKs received from d

# Locks for avoiding race conditions
buffer_lock = Lock()
ack_lock = Lock()


timeout = 0.05
sending_time = 0


# Compute checksum value to append to header
# Returns 32 bytes hash
def checksum(message):
	m = hashlib.md5(message)
	return m.hexdigest()


# Checks if the checksum value in the header matches the checksum value
# computed by the b node.
# Returns true if the received packet is corrupted.
def isCorrupted(message, checksum_header):
	if str(checksum(message)) != str(checksum_header):
		return True
	return False


# Sends the messages received from s to the next hop address.
# Exploits both routers alternatingly (e.g. if previous packet is sent to r1, sent this packet to r2)
def sendMsgToRouters(sockbr1, sockbr2):
	global last_sent_packet, base, window_size, sending_time

	while True:
		buffer_lock.acquire()
		while last_sent_packet < len(send_buffer):				# While there is data to send
			if (base + window_size) > last_sent_packet:			# If number of nonACKed packets is smaller than window size, continue sending
				sending_time = time.time()

				if last_sent_packet % 2 == 1:           		# If last packet is sent to r1
					# Directly send this packet to the next hop(r2)
					sockbr2.sendto(packetizeMsg(last_sent_packet, send_buffer[last_sent_packet]), (ip_r2, port))
					#print "sent to r2"
					print "Packet", last_sent_packet, "is sent."     # For debug
				else:                                   		# If last packet is sent to r2
					# Directly send this packet to the next hop(r1)
					sockbr1.sendto(packetizeMsg(last_sent_packet, send_buffer[last_sent_packet]), (ip_r1,port))
					#print "sent to r1"
					print "Packet", last_sent_packet, "is sent."	 # For debug

				last_sent_packet += 1

			else:												# Wait for timeout
				now = time.time()
				time_passed = now - sending_time				# Calculate time passed after sending
				if time_passed < timeout:						# If it did not wait enough, wait until timeout
					time.sleep(timeout - time_passed)
					#print "sleeeeeeeeeeep"

					if (base + window_size) > last_sent_packet:			# If new ACK received, base is updated
						#print base, window_size, last_sent_packet
						continue
					
				# Retransmit the first nonACKed packet in the window
				ack_lock.acquire()
				packet_num = received_acks.index(0)		# Find the first nonACKed packet
				#print "packet_num", packet_num
				ack_lock.release()
				sockbr2.sendto(packetizeMsg(packet_num, send_buffer[packet_num]), (ip_r2, port))
				sending_time = time.time()
				#print "retransmitted"

		buffer_lock.release()


# Receives messages from s by using TCP socket.
# Target function for receiving thread.
def receiveMsgFromS(client):
	global packet_number

	while True:
		received_message = client.recv(chunk_size, socket.MSG_WAITALL)		# received_message contains the data sent from s

		if received_message:
			print "Packet", packet_number, "is received."

			buffer_lock.acquire()
			send_buffer.append(received_message)							# Store received_message in send_buffer
			buffer_lock.release()

			packet_number += 1



# Receives messages from routers.
# Target function for receiving threads.
def receiveMsgFromRouters(sock):
	global base

	while True:     
		
		received_message, address = sock.recvfrom(chunk_size)				# received_message stores the data received
																			# address stores the (ip, port) tuple of the sender(r1 or r2)
		if received_message:

			check = received_message[0:32]									# Checksum value in header
			pack_num = int(received_message[32:36])							# Packet number in header

			if not isCorrupted(received_message[32:], check):				# Check if ACK packet is corrupted
				ack_lock.acquire()
				received_acks[pack_num] = 1									# If not set ACK is received
			
				base = received_acks.index(0) - 1							# Base is firt nonACKed packet - 1
				ack_lock.release()

			print "ACK", pack_num, "is received:", received_message



# Add sequence number(packet number in our case) to the RDT header
# Add checksum of the packet to RDT header so that 
# when segment is unpacked errors could be detected.       
def packetizeMsg(packet_number, message):
	packnum_str = str(packet_number)
	while len(packnum_str) < 4:				# Since chunks are at most 600 bytes, packet number can have at most 4 digits
											# If less, pad it with zeros
		packnum_str = "0" + packnum_str

	packet = packnum_str + message 			# Add sequence number to the beginning of the message
	
	chk = str(checksum(packet))

	packet = chk + packet 				# Add checksum value to the beginning of the message

	return packet


# Creates sockets for sending and receiving messages.
# Starts the threads for sending and receiving messages.
def connection():
	socksb = socket.socket(socket.AF_INET, socket.SOCK_STREAM)            # Create TCP socket for receiving the message from s
	socksb.bind((ip_s, port))                                             # Bind to s

	sockbr1 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)          # Create UDP socket for sending the message to r1
	sockbr2 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)          # Create UDP socket for sending the message to r2

	sockr1b = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)          # Create UDP socket for receiving the message from r1
	sockr1b.bind((ip_r1_rcv, port)) 

	sockr2b = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)          # Create UDP socket for receiving the message from r2
	sockr2b.bind((ip_r2_rcv, port)) 

	socksb.listen(4)

	client, address = socksb.accept()                                     # Accept the TCP connection request from s

	receive_thread = Thread(target = receiveMsgFromS, args = (client,))	# Create and start thread for receiving from s
	receive_thread.start()
	
	send_thread = Thread(target = sendMsgToRouters, args = (sockbr1, sockbr2))	# Create and start thread for sending to routers
	send_thread.start()

	r1_receive = Thread(target = receiveMsgFromRouters, args = (sockr1b,))	# Create and start thread for receiving from r1
	r1_receive.start()

	r2_receive = Thread(target = receiveMsgFromRouters, args = (sockr2b,))	# Create and start thread for receiving from r2
	r2_receive.start()


if __name__ == "__main__":
	connection()
