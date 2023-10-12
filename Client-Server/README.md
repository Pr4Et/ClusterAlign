# ClusterAlign Client-Server
Alignment of S/TEM tilt series based on clusters of fiducial markers using a client-server approach <br />
Written by Shahar Seifer, in Elbaum lab, Weizmann Institute of Science <br />
Citation: Seifer,S.& Elbaum, M., ClusterAlign: A fiducial tracking and tilt series alignment tool for thick sample tomography. Biological Imaging, 2, E7 (2022). doi:10.1017/S2633903X22000071  <br/>

The client side is run from a single Python module main.py, where the settings are configured in the user manageable sections. It is easy to modify the code for scripts in batch processing. <br />
The client program requires a shared directory, one of the following options:<br />
1. User home folder on the server set by IT to be available within docker containers.<br />
2. Any server that can be accessed by rclone command. You should provide the configuration file set by rclone conf on your own computer.<br />
3. A username and password of FTP site.<br />

Installation of the docker container on the server (or your own computer) is achieved by a single docker command:<br />
docker pull 0525214954/clusteralign_server:latest<br />
Running the docker container on the server:<br />
docker run -it --rm -p <server IP>:110:110 0525214954/clusteralign_server <br />

About<br />
The docker image is adapated from ClusterAlign code, where EMGU-CV was replaced by OPENCVSHARP (https://github.com/shimat/opencvsharp), and installed with all prerequests in the docker image.<br />
The Matlab environement was replaced by Octave that was installed in the docker image with all needed fucntions and Astra Toolbox library.<br />
When running the docker image the server program runs automatically. However, you can stop the program by ctrl+C, edit the source code using nano, complie using dotnet build -f net6.0, and run again using dotnet run.<br />
Changes in the server program may be permanent by ommiting the --rm switch, applying the changes in code, typing exit, and using docker commit.<br/>
The client and server are easy to install in any platform and do not require knowledge in programming.<br />


