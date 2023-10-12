#ClusterAlign client  (ver. 2023Sep30)
#Written by Shahar Seifer (C) 2023, Elbaum lab, Weizmann Institute of Science
#GPL-v3 license
#Tested on PyCharm compiler
#Instructions: edit all the sections in the code written "USER MANAGEABLE".
#For docker installation on server, see https://docs.docker.com/desktop/install/ubuntu/
#The server program is found in docker.io/0525214954/clusteralign_server
#For server installation use: docker pull 0525214954/clusteralign_server:latest
#The file directory may be either the user home directory on the server (ServerLocalDir) or any remote storage configured by rclone.
#Use rclone config to configure a remote file system (example on AWS:  S3:<name of bucket>), then you should write below the path to the rclone.conf file that
# contains your credentials in order to share your remote storage with your own the Docker container.
#For running your docker container on the server use:
# docker run -it --rm -p <server IP>:110:110 0525214954/clusteralign_server
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
import win32api

# Example of setting a server on EC2 AWS (Amazon):
# Create an EC2 instance with sufficient resources (at least 40GB, and several CPU cores), with ports 21 and 110 open in security configurations,
# with Amazon Linux, and click connect
# In Amazon Linux command line: sudo yum install docker,  sudo systemctl start docker,
# docker pull 0525214954/clusteralign_server:latest
# nano .bashrc,   add the line:
# docker run -it --rm -p 110:110 --privileged 0525214954/clusteralign_server
# Stop the instance (to avoid charges).
# Create Elastic IP (permanent address) and attach it to the instance.

context = zmq.Context()
# Socket to talk to server
socket = context.socket(zmq.REQ)
# <<<<<< USER MANAGEABLE :
socket.connect("tcp://132.77.57.166:110")  # write tcp://<IP of server>:110
useFileStorageType=2   #1- local folder on server,  2-remote folder set by rclone, 3-FTP server

if useFileStorageType==3:
    # FTP connection parameters
    ftpHost = '10.0.0.17' # Write IP of FTP server here or DNS name of your organization FTP
    ftpUname = 'ftpuser'
    ftpPass = 'stem'
    ftpPath = '/clusteralign' # actual folder in server is /home/ftpuser' but we write '' if we are using default vsftpd settings
elif useFileStorageType==2:
    RemoteStorage='S3:elbaumlab'  #The name of remote server defined by rclone config : bucket
    LocationRcloneConf=os.getenv('APPDATA')+'\\rclone\\rclone.conf' #Default for windows
    #LocationRcloneConf ='~/.config/rclone.conf' #Default for Linux, Ubuntu
elif useFileStorageType == 1:
    ServerLocalDir='/home/data'
# >>>>>>>>>>>>>>>>>>>>>>


def main():
    #<<<<<< USER MANAGEABLE :
    DataFileName="DM126_2_LT_Cell_1_3.mrc"  #Tag names=file, use filea,fileb for dual axis
    TiltFileName="DM126_2_LT_Cell_1_3.rawtlt"
    if useFileStorageType==3:
        LocalPath = "e:/volume"  # your local path where the data is found
        LocalDataFileName = "DM126_2_LT_Cell_1_3.mrc"  # Your local dataset filename
        LocalTiltFileName = "DM126_2_LT_Cell_1_3.rawtlt"  # your local angle filename
        # >>>>>>>>>>>>>>>>>>>>>>
        showFTP(ftpHost,ftpUname,ftpPass,ftpPath)
        purgeFTP(ftpHost,ftpUname,ftpPass,ftpPath)  #Remove previous files from FTP folder
        sendfile(LocalPath,LocalDataFileName,DataFileName,ftpHost,ftpUname,ftpPass,ftpPath)  #send to FTP
        sendfile(LocalPath,LocalTiltFileName,TiltFileName,ftpHost,ftpUname,ftpPass,ftpPath)

    if useFileStorageType == 2:
        say("Rclone")
        #Send file content of LocationRcloneConf, then send the word Endconfig
        file1 = open(LocationRcloneConf, 'r')
        Lines = file1.readlines()
        for line in Lines:
            say(line)
        say("Endconfig")
    #Take Server
    say("Set") #Mark the start of alignment configuration
    if useFileStorageType == 3:
        say("ftpPath={}".format(ftpPath))
        say("ftpHost={}".format(ftpHost))
        say("ftpUname={}".format(ftpUname))
        say("ftpPass={}".format(ftpPass))
        say("useFTP=1")
    elif useFileStorageType == 2:
        say("RemoteMount={}".format(RemoteStorage))
    elif useFileStorageType == 1:
        say("ServerLocalDir={}".format(ServerLocalDir))


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

    if useFileStorageType == 3:
        getfile(LocalPath, os.path.splitext(DataFileName)[0] + ".output.txt",
                os.path.splitext(DataFileName)[0] + ".output.txt", ftpHost, ftpUname, ftpPass, ftpPath)
        f = open(LocalPath + "/" + os.path.splitext(LocalDataFileName)[0] + ".output.txt", "r")
        print(f.read())
        f.close()
        #Download the final alignement from the FTP to your local folder
        getfile(LocalPath, os.path.splitext(LocalDataFileName)[0]+"_ali.mrc", os.path.splitext(DataFileName)[0]+"_ali.mrc", ftpHost, ftpUname, ftpPass, ftpPath)
    elif useFileStorageType == 1:
        f = open(ServerLocalDir + "/" + os.path.splitext(LocalDataFileName)[0] + ".output.txt", "r")
        print(f.read())
        f.close()
        print("The output files are on your home directory")
    elif useFileStorageType == 2:
        print("The output files are on the file server")

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


if __name__=="__main__":
    main()


