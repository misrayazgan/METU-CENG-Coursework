ip = "10.10.1.2"            # ip address of the next hop
port = 12346
measurements = ["Current temperature: 23 C", "Current temperature: 24 C", "Current temperature: 25 C",
"Current temperature: 26 C", "Current temperature: 27 C", "Current temperature: 28 C",
"Current temperature: 29 C", "Current temperature: 30 C", "Current temperature: 31 C",
"Current temperature: 32 C"]

messages = measurements * 10

chunk_size = 39			# To send the messages in chunks of 39 bytes

import socket
from datetime import datetime


# Add timestamps to the messages to indicate the sending time.
# Timestamps will be used later in the d node to calculate the end-to-end delay
def packetizeMsg(message):
    cur_time = datetime.now()
    msec_time = int(float(cur_time.strftime("%s.%f")) * 1000)
    return message + ';' + str(msec_time)						# Use ';' as a delimiter to distinct the measurement result from the timestamp


# Creates the socket and sends the message to the next hop address
def connection():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)    # Create tcp socket
    sock.connect((ip, port))                                    # Connect the socket to b's address (start the three-way handshake)

    # Send all the sensor measurements to b
    packet_number = 1
    for message in messages:
        sock.send(packetizeMsg(message))
        print "Packet", packet_number, "is sent, the content is: " packetizeMsg(message)		# For debug
        packet_number += 1



if __name__ == "__main__":
    connection()
