/*
 * Name: Misranur
 * Surname: Yazgan
 *
 */

Firstly, we placed each of the scripts to the related resource nodes which are reserved on the Geni Portal. From 5 different terminal windows we connected to 
each of the nodes by typing the following command: ssh -i <private key location> <username>@<hostname> -p <portnumber> 
where port numbers are 28710, 28706, 28708, 28709 and 28707 for s, b, r1, r2 and d respectively.


***HOW TO RUN:

There are 5 scripts named s.py, b.py, r1.py, r2.py and d.py, each running on the nodes s, b, r1, r2 and d correspondingly.
The script for the Broker node(b.py) is executed firstly as it is stated in the homework text, then r1.py, r2.py, d.py and s.py are executed respectively.
Only s.py completes its execution after sending the messages whereas b.py, r1.py, r2.py and d.py continue to listen and do not finish their executions.

In order to stop their executions manually, keyboard interrupt (ctrl + C) can be used.

***Message format: The messages to be sent are given with a list inside s.py. We conducted our experiments with 100 packets(messages) each of which are 39 bytes long and 
in the "Current temperature: <integer> C;<timestamp>" format.


***To synchronize the nodes we used the following commands in both s and d nodes:

sudo service ntp stop
sudo ntpdate -s time.nist.gov
sudo service ntp start



***Netem/tc Link Configuration commands:

For 1ms +-5ms network emulation delay:
1. On the b node: eth2 is for the link between b and r1, eth3 is for the link between b and r2
	sudo tc qdisc add dev eth2 root netem delay 1ms 5ms distribution normal		(since this is the first command for the delay on link eth2, we used add)
	sudo tc qdisc add dev eth3 root netem delay 1ms 5ms distribution normal		(since this is the first command for the delay on link eth3, we used add)

2. On the r1 node: eth2 is for the link between r1 and d
	sudo tc qdisc add dev eth2 root netem delay 1ms 5ms distribution normal		(since this is the first command for the delay on link eth2, we used add)

3. On the r2 node: eth2 is for the link between r2 and d
	sudo tc qdisc add dev eth2 root netem delay 1ms 5ms distribution normal		(since this is the first command for the delay on link eth2, we used add)


For 20ms +-5ms network emulation delay:
1. On the b node:
	sudo tc qdisc change dev eth2 root netem delay 20ms 5ms distribution normal		(since we used add for eth2 in the 1ms part, we can use change)
	sudo tc qdisc change dev eth3 root netem delay 20ms 5ms distribution normal		(since we used add for eth3 in the 1ms part, we can use change)

2. On the r1 node:
	sudo tc qdisc change dev eth2 root netem delay 20ms 5ms distribution normal		(since we used add for eth2 in the 1ms part, we can use change)

3. On the r2 node:
	sudo tc qdisc change dev eth2 root netem delay 20ms 5ms distribution normal		(since we used add for eth2 in the 1ms part, we can use change)


For 60ms +-5ms network emulation delay:
1. On the b node:
	sudo tc qdisc change dev eth2 root netem delay 60ms 5ms distribution normal		(since we used add for eth2 in the 1ms part, we can use change)
	sudo tc qdisc change dev eth3 root netem delay 60ms 5ms distribution normal		(since we used add for eth3 in the 1ms part, we can use change)

2. On the r1 node:
	sudo tc qdisc change dev eth2 root netem delay 60ms 5ms distribution normal		(since we used add for eth2 in the 1ms part, we can use change)

3. On the r2 node:
	sudo tc qdisc change dev eth2 root netem delay 60ms 5ms distribution normal		(since we used add for eth2 in the 1ms part, we can use change)



***NOTE: We are providing the figures included in the tex file for convenience.
