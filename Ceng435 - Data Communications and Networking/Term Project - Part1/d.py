ip_r1 = "10.10.3.2"            # ip address of the d-interface for receiving a message from r1
ip_r2 = "10.10.5.2"            # ip address of the d-interface for receiving a message from r2
port = 12346
chunk_size = 39				   # To send the messages in chunks of 39 bytes

import socket
import select
from datetime import datetime


# Unpacks the packet according to the message format and gets the timestamp of the sending time from s.
# Then, calculates the end-to-end delay and prints it to use it in the experiment.
def totalTranmissionTime(received_message, packet_number):
    cur_time = datetime.now()
    msec_time = int(float(cur_time.strftime("%s.%f")) * 1000)
    print packet_number, msec_time - int(received_message[26:])


# Receives messages from r1 and r2, calls totalTranmissionTime() to calculate the delay
def receiveMsg(sockets):
    packet_number = 1
    while True:
        read_ready, write_ready, exc = select.select(sockets, [], [])       # Select from which socket to read

        for s in read_ready:
            received_message, address = s.recvfrom(chunk_size)              # received_message stores the data received
                                                                    		# address stores the (ip, port) tuple of the sender(r2)
            if received_message:
                totalTranmissionTime(received_message, packet_number)
            packet_number += 1


# Create UDP socket for receiving messages from the routers (r1 or r2)
def connection():
    sockr1d = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)             # Create a UDP socket for receiving messages from r1
    sockr1d.bind((ip_r1, port))                                            # Bind the socket to r1's address

    sockr2d = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)             # Create a UDP socket for receiving messages from r2
    sockr2d.bind((ip_r2, port))                                            # Bind the socket to r2's address

    sockets = [sockr1d, sockr2d]

    receiveMsg(sockets)



if __name__ == "__main__":
    connection()
