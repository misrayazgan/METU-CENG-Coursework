#include <stdio.h>
#include <pcap.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#define SIZE_ETHERNET 14


// Struct for IP header
struct ipHeader
{
    u_char  ip_vhl;                 // IP version
    u_char  ip_tos;                 // Type of service
    u_short ip_len;                 // Total length
    u_short ip_id;                  // Identification
    u_short ip_off;                 // Fragment offset field
    #define IP_RF 0x8000            // Reserved fragment flag
    #define IP_DF 0x4000            // Dont fragment flag
    #define IP_MF 0x2000            // More fragments flag
    #define IP_OFFMASK 0x1fff       // Mask for fragmenting bits
    u_char  ip_ttl;                 // Time to live
    u_char  ip_p;                   // Protocol
    u_short ip_sum;                 // Checksum
    struct  in_addr ip_src,ip_dst;  // Source and dest address
};


// Struct for TCP header
struct tcpHeader
{
    u_short th_sport;               // Source port
    u_short th_dport;               // Destination port
    u_int th_seq;                   // Sequence number
    u_int th_ack;                   // Acknowledgement number
    u_char  th_offx2;               // Data offset, rsvd
	#define TH_OFF(th)      (((th)->th_offx2 & 0xf0) >> 4)
    u_char  th_flags;
    #define TH_FIN  0x01
    #define TH_SYN  0x02
    #define TH_RST  0x04
    #define TH_PUSH 0x08
    #define TH_ACK  0x10
    #define TH_URG  0x20
    #define TH_ECE  0x40
    #define TH_CWR  0x80
    #define TH_FLAGS        (TH_FIN|TH_SYN|TH_RST|TH_ACK|TH_URG|TH_ECE|TH_CWR)
    u_short th_win;                 // Window
    u_short th_sum;                 // Checksum
    u_short th_urp;                 // Urgent pointer
};



// Helper function for printing the lines of the payload.
void printLineContent(const u_char *payload, int length, int offset)
{
	int i;

	// If valid ASCII character, print the content of payload
	for(i = 0; i < length; i++)
	{
		if(*payload >= 0 && *payload <= 127)
			printf("%c ", *payload);
		payload++;
	}

	printf("\n");

	return;
}


// Function for printing the payload of a packet.
void printPayload(const u_char *payload, int length)
{
	int lineWidth = 16;
	int lineLength;
	int offset = 0;

	if(length <= 0)
		return;


	if(length <= lineWidth)
	{
		printLineContent(payload, length, offset);
		return;
	}


	while(1)
	{
		lineLength = lineWidth % length;

		printLineContent(payload, lineLength, offset);

		length -= lineLength;

		payload += lineLength;			// shift pointer

		offset += lineWidth;

		if(length < lineWidth)
		{
			printLineContent(payload, length, offset);
			break;
		}
	}

	return;
}


// Function for handling captured packets.
void packetCaptured(u_char *args, const struct pcap_pkthdr *header, const u_char *packet)
{
	const struct ipHeader *ip;				// IP header
	const struct tcpHeader *tcp;			// TCP header
	static int packetNumber = 1;

	int sizeTCP;
	int sizeIP;
	
	// Print number of captured packet
	printf("\nPacket number: %d\n", packetNumber);
	packetNumber++;

	ip = (struct ipHeader*)(packet + SIZE_ETHERNET);
	sizeIP = ((ip->ip_vhl) & 0x0f) * 4;

	if(sizeIP < 20)
	{
		printf("Invalid IP header length: %u bytes\n", sizeIP);
		return;
	}


	/*****Print source and destination IP addresses*****/
	printf("Source IP: %s\n", inet_ntoa(ip->ip_src));
	printf("Destination IP: %s\n", inet_ntoa(ip->ip_dst));


	if(ip->ip_p == IPPROTO_IP)
	{
		printf("Protocol: IP\n");
	}
	else if(ip->ip_p == IPPROTO_ICMP)
	{
		printf("Protocol: ICMP\n");
	}
	else if(ip->ip_p == IPPROTO_TCP)			// Print TCP packet's payload
	{
		printf("Protocol: TCP\n");

		const u_char *payload;
		int sizePayload;

		tcp = (struct tcpHeader*)(packet + SIZE_ETHERNET + sizeIP);
		sizeTCP = TH_OFF(tcp) * 4;

		if(sizeTCP < 20)
		{
			printf("Invalid TCP header length: %u bytes\n", sizeTCP);
			return;
		}


		printf("Source port: %d\n", ntohs(tcp->th_sport));
		printf("Destination port: %d\n", ntohs(tcp->th_dport));


		payload = (u_char*)(packet + SIZE_ETHERNET + sizeIP + sizeTCP);

		sizePayload = ntohs(ip->ip_len) - (sizeIP + sizeTCP);

		if(sizePayload > 0)
		{
			printPayload(payload, sizePayload);
		}
	}

	return;
}



int main()
{
	char *device;
	char errorBuffer[PCAP_ERRBUF_SIZE];
	pcap_t *handle;
	struct bpf_program fp;			// compiled filter expression

	/*******Filter expressions*******/
	char filterExp[] = "tcp dst portrange 10-100";
	//char filterExp[] = "icmp and host 10.0.2.4";
	//char filterExp[] = "icmp and host 10.0.2.5 and host 10.0.2.4";

	bpf_u_int32 subnetMask;			// our subnet mask
	bpf_u_int32 ourIP;				// our IP
	int numberOfPackets = 10000;


	// Look for the device name
	device = pcap_lookupdev(errorBuffer);
	if(device == NULL)
	{
		fprintf(stderr, "Could not find a device: %s\n", errorBuffer);
		exit(EXIT_FAILURE);
	}


	// Get subnet mask and IP of the device
	if(pcap_lookupnet(device, &ourIP, &subnetMask, errorBuffer) == -1)
	{
		fprintf(stderr, "Cannot get subnetMask: %s\n", errorBuffer);
		subnetMask = 0;
		ourIP = 0;
	}


	// Open capture device
	handle = pcap_open_live(device, BUFSIZ, 1, 1000, errorBuffer);
	if(handle == NULL)
	{
		fprintf(stderr, "Could not open device: %s\n", errorBuffer);
		exit(EXIT_FAILURE);
	}


	// Make sure about capturing on an Ethernet device
	if(pcap_datalink(handle) != DLT_EN10MB)
	{
		fprintf(stderr, "Ethernet header is not supported: %s\n", errorBuffer);
		exit(EXIT_FAILURE);
	}


	// Compile the filter expression
	if(pcap_compile(handle, &fp, filterExp, 0, ourIP) == -1)
	{
		fprintf(stderr, "Could not parse filter: %s\n", pcap_geterr(handle));
		exit(EXIT_FAILURE);
	}


	// Set the compiled filter
	if(pcap_setfilter(handle, &fp) == -1)
	{
		fprintf(stderr, "Could not install filter: %s\n", pcap_geterr(handle));
		exit(EXIT_FAILURE);
	}


	// Set callback function
	pcap_loop(handle, numberOfPackets, packetCaptured, NULL);

	// Clean and exit
	pcap_freecode(&fp);
	pcap_close(handle);

	return 0;
}
