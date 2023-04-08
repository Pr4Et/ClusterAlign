#ClusterAlign client  (ver. 2023Apr08)
#Written by Shahar Seifer (C) 2023, Elbaum lab, Weizmann Institute of Science
#GPL-v3 license
#Tested on PyCharm compiler
#For docker installation on server, see https://docs.docker.com/desktop/install/ubuntu/
#For docker installation on AWS EC2 server using Amazon Linux:  sudo yum install docker,  sudo systemctl start docker,  (open ports 21 and 110 in security configurations).
#The server program is found in docker.io/0525214954/clusteralign_server
#For server installation use: docker pull 0525214954/clusteralign_server:latest
#For setting up an FTP server (typically on the server) see details in the end of this code, or use: https://www.geeksforgeeks.org/how-to-setup-and-configure-an-ftp-server-in-linux-2/ avoiding security measures
#For running docker container on the server use:
# docker run -it --rm -p <server IP>:110:110 --privileged 0525214954/clusteralign_server
#To access the server source code or run Octave scripts you can stop the docker container by Ctrl+C on the server side.
#Article/citation: Seifer, S., & Elbaum, M. (2022). ClusterAlign: A fiducial tracking and tilt series alignment tool for thick sample tomography. Biological Imaging, 2, E7. doi:10.1017/S2633903X22000071

import os
import string
from PIL import Image
import zmq
import ftplib
import base64
import io
import psutil
#When using EC2 server on AWS (alternative to local server) you need the following to start and stop the container from running:
#import sys
#import boto3
#from botocore.exceptions import ClientError
#ec2 = boto3.client('ec2')
#EC2_instance_id=


context = zmq.Context()
# Socket to talk to server
socket = context.socket(zmq.REQ)
socket.connect("tcp://10.0.0.17:110")  # write tcp://<IP of server>:110
# FTP connection parameters
ftpHost = '10.0.0.17' # Write IP of FTP server here or DNS name of your organization FTP
ftpUname = 'ftpuser'
ftpPass = 'stem'
ftpPath = '' # actual folder in server is /home/ftpuser' but we write '' if we are using default vsftpd settings


def main():
    #<<<<<< USER MANAGEABLE :
    operation_step = 0  #0= alignment, 1= reconstruction, 2= cleaning and stopping
    LocalPath="e:/volume" #your local path where the data is found
    LocalDataFileName="DM126_2_LT_Cell_1_3.mrc" # Your local dataset filename
    LocalTiltFileName="DM126_2_LT_Cell_1_3.rawtlt"  #your local angle filename
    # >>>>>>>>>>>>>>>>>>>>>>

    DataFileName="file.mrc"  #Tag names=file, use filea,fileb for dual axis
    TiltFileName="file.rawtlt"
    showFTP(ftpHost,ftpUname,ftpPass,ftpPath)
    #Alignment
    if operation_step == 0:
        #AWS_EC2(EC2_instance_id, "ON")  #Start to run container on AWS server (charges may apply). Alternative for local server
        purgeFTP(ftpHost,ftpUname,ftpPass,ftpPath)  #Remove previous files from FTP folder
        sendfile(LocalPath,LocalDataFileName,DataFileName,ftpHost,ftpUname,ftpPass,ftpPath)  #send to FTP
        sendfile(LocalPath,LocalTiltFileName,TiltFileName,ftpHost,ftpUname,ftpPass,ftpPath)

        #Take Server
        say("Set") #Mark the start of alignment configuration
        say("ftpPath={}".format(ftpPath))
        say("ftpHost={}".format(ftpHost))
        say("ftpUname={}".format(ftpUname))
        say("ftpPass={}".format(ftpPass))

        # <<<<<< USER MANAGEABLE :
        #alignment parameters
        say("DataFileName={}".format(DataFileName)) #Name of mrc file in FTP site, without path
        say("TiltFileName={}".format(TiltFileName)) #Name of rawtlt file in FTP site, without path
        say("fiducials_bright=False") #Dark fiducials should be marked with False, otherwise =True.
        say("xisRotation=True") #=True for horizontal axis of rotation, =False for vertical, as you observe approximatley from the stack
        say("cluster_size=400") #Radius in pixels for collecting vectors between fiducial particles
        say("NfidMax=800") #Maximum number of particles to look for. Specify about twice the expected number of particles in the image to have low false negative
        say("fidsize=7.5") #Average size of particles, in pixels
        say("minimum_tracked_fiducials=90") #number in percents: filter fiducials that do not occur in as match percent of tilt views
        say("PreAlignmentTol=200") #You set the maximum shift allowed for alignement (in pixels)
        say("TolFidCenter=200") #Set between 200 and 300 (in percents compared to fiducial size, the allowed distortion in affinity of vectors between particles)
        say("TolFidSize=30") #Variability (in percents) in fiducial particles sizes
        say("Ncluster=-1") #Leave "Ncluster=-1", or specify number of vectors necessary to compare cluster signatures
        say("ncenter=-1") #Leave "ncenter=-1" for automatic finding the center slice, or enter the slice number of the zero tilt
        say("isArbAngle=False") #Usually we set "isArbAngle=False", unless the rotation axis is not close to 0 or 90 degrees
        say("coswindow=False") #Usually you set "coswindow=False", unless the scan is based on Hoppe's cos samping grid
        say("export_normalali=True") #Set "export_normalali=True" for generating ali files compatible with IMOD
        #say("fidfileName=") #Either leave "fidfileName=" or specify previous .fid file to start from
        # >>>>>>>>>>>>>>>>>>>>>>

        say("Start")
        print("Alignment started")
        message=""
        while (message.startswith("see") or message=="" or message.startswith("calc")):
            socket.send_string("OK", 0)
            message = socket.recv_string(0)
            if message.startswith("see"):
                base64_bytes = message[3:].encode('ascii')
                imgdata=base64.decodebytes(base64_bytes)
                ImID=Image.open(io.BytesIO(imgdata)).show()
            else:
                print(message)

        getfile(LocalPath, os.path.splitext(LocalDataFileName)[0]+".output.txt", os.path.splitext(DataFileName)[0]+".output.txt", ftpHost, ftpUname, ftpPass, ftpPath)
        f = open(LocalPath+"/" +os.path.splitext(LocalDataFileName)[0]+".output.txt", "r")
        print(f.read())
        f.close()
        #Download the final alignement from the FTP to your local folder
        getfile(LocalPath, os.path.splitext(LocalDataFileName)[0]+"_ali.mrc", os.path.splitext(DataFileName)[0]+"_ali.mrc", ftpHost, ftpUname, ftpPass, ftpPath)

    #Reconstruction
    if operation_step == 1:
        say("Set") #Mark the start of configuration
        say("DataFileName={}".format(DataFileName))  # Name of mrc file in FTP site, without path

        # <<<<<< USER MANAGEABLE :
        say("Reconstruct_bin=4")
        say("Reconstruct_thickness=250")
        #>>>>>>>>>>>>>>>>>>>>>>>>

        say("Reconstruct")
        socket.send_string("OK", 0)
        message = socket.recv_string()
        print(message)

    #Cleaning and stopping
    if operation_step == 2:
        say("Purge")  #remove files from local folder (/shared) on server
        # AWS_EC2(EC2_instance_id, "OFF")  #Stop container on AWS server (you can also set automatic termination in their site)

    socket.close()


#send your local files to your FTP server
def sendfile(LocalPath,Localfilename,filename,ftpHost,ftpUname,ftpPass,ftpPath):
    print(f'Sending {Localfilename} to FTP')
    ftpPort = 21
    localFilePath = f'{LocalPath}/{Localfilename}'
    ftp = ftplib.FTP(timeout=12)
    ftp.connect(ftpHost, ftpPort)
    # login to the FTP server
    ftp.login(ftpUname, ftpPass)
    ftp.cwd(ftpPath)
    # Read file in binary mode
    with open(localFilePath, "rb") as file:
        # upload file to FTP server using storbinary, specify blocksize(bytes) only if higher upload chunksize is required
        retCode = ftp.storbinary(f"STOR {filename}", file, blocksize=1024 * 1024)
    ftp.quit()

def getfile(LocalPath,Localfilename,filename,ftpHost,ftpUname,ftpPass,ftpPath):
    print(f'Receiving {Localfilename} from FTP')
    ftpPort = 21
    localFilePath = f'{LocalPath}/{Localfilename}'
    ftp = ftplib.FTP(timeout=12)
    ftp.connect(ftpHost, ftpPort)
    # login to the FTP server
    ftp.login(ftpUname, ftpPass)
    ftp.cwd(ftpPath)
    # Read file in binary mode
    with open(localFilePath, "wb") as file:
        # use FTP's RETR command to download the file
        ftp.retrbinary(f"RETR {filename}", file.write, blocksize=1024 * 1024)
    ftp.quit()

def purgeFTP(ftpHost,ftpUname,ftpPass,ftpPath):
    print('Deleting files on FTP folder')
    ftpPort = 21
    ftp = ftplib.FTP(timeout=12)
    ftp.connect(ftpHost, ftpPort)
    # login to the FTP server
    ftp.login(ftpUname, ftpPass)
    ftp.cwd(ftpPath)
    for something in ftp.nlst():
        try:
            ftp.delete(something)
        except Exception:
            print('')
    ftp.quit()

def showFTP(ftpHost,ftpUname,ftpPass,ftpPath):
    print('Files on FTP folder:')
    ftpPort = 21
    ftp = ftplib.FTP(timeout=12)
    ftp.connect(ftpHost, ftpPort)
    # login to the FTP server
    ftp.login(ftpUname, ftpPass)
    ftp.cwd(ftpPath)
    for something in ftp.nlst():
        print(something)
    ftp.quit()



def say(what):
    socket.send_string(what,0)#flags: 0, NOBLOCK, SNDMORE, or NOBLOCK|SNDMORE
    #  Get the reply.
    message = socket.recv_string(0) #flags (int) – 0 or NOBLOCK.
    print(message)

def AWS_EC2(instance_id,action):
    #Copied from https://boto3.amazonaws.com/v1/documentation/api/latest/guide/ec2-example-managing-instances.html
    if action == 'ON':
        # Do a dryrun first to verify permissions
        try:
            ec2.start_instances(InstanceIds=[instance_id], DryRun=True)
        except ClientError as e:
            if 'DryRunOperation' not in str(e):
                raise

        # Dry run succeeded, run start_instances without dryrun
        try:
            response = ec2.start_instances(InstanceIds=[instance_id], DryRun=False)
            print(response)
        except ClientError as e:
            print(e)
    else:
        # Do a dryrun first to verify permissions
        try:
            ec2.stop_instances(InstanceIds=[instance_id], DryRun=True)
        except ClientError as e:
            if 'DryRunOperation' not in str(e):
                raise

        # Dry run succeeded, call stop_instances without dryrun
        try:
            response = ec2.stop_instances(InstanceIds=[instance_id], DryRun=False)
            print(response)
        except ClientError as e:
            print(e)


if __name__=="__main__":
    main()


"""
How to set up FTP server with username and password that both server and client can access
Important notes:
1. use dedicated folder for the clusteralign files, since purge command will remove all files on the directory
2. In case your organization blocks file sharing you are encouraged to use your orginzation dedicated FTP site

Instructions for FTP server settings on Ubuntu: 

sudo apt install vsftpd

(To verify status is active:)
sudo systemctl status vsftpd

(To allow in Firewall:)
sudo ufw allow 20/tcp
sudo ufw allow 21/tcp

sudo adduser ftpuser
(enter any password for user ftpuser and skip the user profile questions)

sudo chown <Linux user> /home/ftpuser

sudo nano /etc/vsftpd.conf
(change to write_enable=YES,  ^o, enter, ^x)
sudo systemctl restart --now vsftpd

To manage the ftp site manually we can type:
ftp ftp://<server IP address>
(type the ftpuser password afterword)
(you can browse the files by ls, add files or download using put and get commands, quit by bye command. see help for other options)

In this Python program fill out the FTP username: ftpuser, and the password you entered, and the IP address/ dns name of the computer.


"""