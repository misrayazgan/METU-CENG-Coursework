/*
 * Name: Misranur
 * Surname: Yazgan
 *
 */



Firstly, we placed each of the scripts to the related resource nodes which are reserved on the Geni Portal. From 5 different terminal windows we connected to each of the nodes by typing the following command: 

ssh -i <private key location> <username>@<hostname> -p <portnumber> 

with the corresponding port numbers for each node.


*****HOW TO RUN*****
There are 3 scripts, s.py in Source node, b.py in Broker node and d.py in Destination node. After conneting to GENI nodes, 
we first run b.py, then d.py and s.py respectively. s.py finishes its execution after sending all the chunks, but b.py and d.py continue listening.

In order to stop their executions manually, keyboard interrupt (ctrl + Z) can be used.


***NOTE: We are also providing the figures that are included in the tex file in zip file for convenience while reading the report.


*** Chunk size: We predefined our chunk size as 600 bytes. After adding header for RDT, it became:
checksum(32 bytes) + packet number(4 bytes) + original message(600 bytes)

So, with header the total packet size becomes 636 bytes which is less than max packet size limit.


*** To synchronize the nodes we used the following commands in both s and d nodes:

sudo service ntp stop
sudo ntpdate -s time.nist.gov
sudo service ntp start


*** Routing Commands:
In node r1, we used the following commands:

sudo route add -net 10.10.3.2 netmask 255.255.255.255 gw 10.10.3.2		(to forward packets coming from b to d)
sudo route add -net 10.10.2.1 netmask 255.255.255.255 gw 10.10.2.1		(to forward packets coming from d to b)

In node r2, we used the following commands:

sudo route add -net 10.10.5.2 netmask 255.255.255.255 gw 10.10.5.2		(to forward packets coming from b to d)
sudo route add -net 10.10.4.1 netmask 255.255.255.255 gw 10.10.4.1		(to forward packets coming from d to b)


*** Routing Table In Node R1:
Destination     Gateway         Genmask         Flags Metric Ref    Use Iface
0.0.0.0         172.16.0.1      0.0.0.0         UG    0      0        0 eth0
10.0.0.0        10.10.2.1       255.0.0.0       UG    0      0        0 eth2
10.10.2.0       0.0.0.0         255.255.255.0   U     0      0        0 eth2
10.10.2.1       10.10.2.1       255.255.255.255 UGH   0      0        0 eth2
10.10.3.0       0.0.0.0         255.255.255.0   U     0      0        0 eth1
10.10.3.2       10.10.3.2       255.255.255.255 UGH   0      0        0 eth1
10.10.4.0       10.10.2.1       255.255.255.0   UG    0      0        0 eth2
10.10.4.0       10.10.3.2       255.255.252.0   UG    0      0        0 eth1
10.10.4.2       10.10.3.2       255.255.255.254 UG    0      0        0 eth1
172.16.0.0      0.0.0.0         255.240.0.0     U     0      0        0 eth0


*** Routing Table in Node R2:
Destination     Gateway         Genmask         Flags Metric Ref    Use Iface
0.0.0.0         172.16.0.1      0.0.0.0         UG    0      0        0 eth0
10.0.0.0        10.10.4.1       255.0.0.0       UG    0      0        0 eth1
10.10.3.2       10.10.5.2       255.255.255.254 UG    0      0        0 eth2
10.10.4.0       0.0.0.0         255.255.255.0   U     0      0        0 eth1
10.10.4.1       10.10.4.1       255.255.255.255 UGH   0      0        0 eth1
10.10.5.0       0.0.0.0         255.255.255.0   U     0      0        0 eth2
10.10.5.2       10.10.5.2       255.255.255.255 UGH   0      0        0 eth2
172.16.0.0      0.0.0.0         255.240.0.0     U     0      0        0 eth0




***Netem/tc Link Configuration commands:

For packet loss experiment 0.5%:
1. On b node: eth2 is for the link between b and r1, eth3 is for the link between b and r2
	sudo tc qdisc add dev eth2 root netem loss 0.5% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%
	sudo tc qdisc add dev eth3 root netem loss 0.5% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

2. On r1 node: eth1 is for the link between r1 and d
	sudo tc qdisc add dev eth1 root netem loss 0.5% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

3. On r2 node: eth2 is for the link between r2 and d
	sudo tc qdisc add dev eth2 root netem loss 0.5% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

(since these are the first commands for configuring the packet loss on corresponding links, we used add, after this we used change command)


For packet loss experiment 10%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 10% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%
	sudo tc qdisc change dev eth3 root netem loss 10% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 10% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 10% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%


For packet loss experiment 20%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 20% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%
	sudo tc qdisc change dev eth3 root netem loss 20% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 20% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 20% corrupt 0% duplicate 0% delay 3ms reorder 0% 0%


For corruption experiment 0.2%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0.2% duplicate 0% delay 3ms reorder 0% 0%
	sudo tc qdisc change dev eth3 root netem loss 0% corrupt 0.2% duplicate 0% delay 3ms reorder 0% 0%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 0% corrupt 0.2% duplicate 0% delay 3ms reorder 0% 0%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0.2% duplicate 0% delay 3ms reorder 0% 0%


For corruption experiment 10%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 10% duplicate 0% delay 3ms reorder 0% 0%
	sudo tc qdisc change dev eth3 root netem loss 0% corrupt 10% duplicate 0% delay 3ms reorder 0% 0%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 0% corrupt 10% duplicate 0% delay 3ms reorder 0% 0%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 10% duplicate 0% delay 3ms reorder 0% 0%


For corruption experiment 20%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 20% duplicate 0% delay 3ms reorder 0% 0%
	sudo tc qdisc change dev eth3 root netem loss 0% corrupt 20% duplicate 0% delay 3ms reorder 0% 0%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 0% corrupt 20% duplicate 0% delay 3ms reorder 0% 0%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 20% duplicate 0% delay 3ms reorder 0% 0%


For reorder experiment 1% 50%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 1% 50%
	sudo tc qdisc change dev eth3 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 1% 50%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 1% 50%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 1% 50%


For reorder experiment 10% 50%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 10% 50%
	sudo tc qdisc change dev eth3 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 10% 50%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 10% 50%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 10% 50%


For reorder experiment 35% 50%:
1. On b node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 35% 50%
	sudo tc qdisc change dev eth3 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 35% 50%

2. On r1 node:
	sudo tc qdisc change dev eth1 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 35% 50%

3. On r2 node:
	sudo tc qdisc change dev eth2 root netem loss 0% corrupt 0% duplicate 0% delay 3ms reorder 35% 50%
