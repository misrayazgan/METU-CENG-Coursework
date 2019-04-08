ip = "10.10.1.2"            # ip address of the b-interface for receiving a message from s
ip_r1 = "10.10.2.2"         # ip address of the next hop(r1)
ip_r2 = "10.10.4.2"         # ip address of the next hop(r2)
port = 12346

chunk_size = 39				# To send the messages in chunks of 39 bytes

import socket
import threading
import sys
from datetime import datetime


# Receives messages from s, sends them to the next hop address (r1 or r2)
# in an alternating way by using the turn variable
def receiveMsg(client, sockbr1, sockbr2, turn):
    received_message = client.recv(chunk_size)                    # received_message contains the data sent from s
    if received_message:
        print "Msg is received, forwarding: ", received_message
        if turn % 2 == 0:
            sockbr2.sendto(received_message, (ip_r2, port))       # Directly send the data to the next hop(r2)
        else:
            sockbr1.sendto(received_message, (ip_r1,port))        # Directly send the data to the next hop(r1)


# Create TCP socket for receiving messages from s
# Create UDP socket to send the messages to the next hop address
def connection():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)            # Create TCP socket for receiving the message from s
    sock.bind((ip, port))                                               # Bind to s

    sockbr1 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)          # Create UDP socket for sending the message to r1
    sockbr2 = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)          # Create UDP socket for sending the message to r2

    sock.listen(4)

    turn = 0                    # When this number is even, send to r2
                                # When odd, send to r1

    while True:
    	client, address = sock.accept()                						# Accept the TCP connection request from s
    	
        receiveMsg(client, sockbr1, sockbr2, turn)
        turn += 1


if __name__ == "__main__":
    connection()
