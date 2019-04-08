ip = "10.10.4.2"            # ip address of the r1-interface for receiving a message from b
ip_d = "10.10.5.2"          # ip address of the next hop
port = 12346

chunk_size = 39				# To send the messages in chunks of 39 bytes

import socket
import sys
from datetime import datetime




def receiveMsg(client, sockr2d):
    received_message, address = client.recvfrom(chunk_size)            # received_message stores the data received
                                                                	   # address stores the (ip, port) tuple of the sender(b)
    print received_message
    if received_message:
        print "Msg is received, forwarding: ", received_message		   # For debug
        sockr2d.sendto(received_message, (ip_d,port))				   # Directly send the message to d



# Create UDP socket for receiving messages from b
# Create UDP socket to send the messages to the next hop address
def connection():
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)            # Create UDP socket for receiving the message from b
    sock.bind((ip, port))                                              # Bind to s

    sockr2d = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)         # Create UDP socket for sending the message to d

    while True:
        receiveMsg(sock, sockr2d)

if __name__ == "__main__":
    connection()
