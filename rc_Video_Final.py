import picamera, os
import datetime as dt
from time import strftime
import subprocess
import time
from subprocess import Popen

#===========================================
# Configuration
#===========================================
# change parameters here before starting program!!

# change to include female IDs: #, type, color e.g. 13VTR
cameraID = "pi_no1"
audioID = ""

female = ""
date=dt.today().strftime('%Y-%m-%d')

#type of observation we are recording
obs = "phonotaxis" 


# video parameters
vidHours = 24 # in hours
segLength = 60 # in minutes
movie_path = '/home/pi/Desktop/Video'


# When to start recording
same_day = True # True (start recording same day) or False (start recording the next day)
start_hour = 18 # Use military time
start_min = 0

print('HELLO')

#===========================================
# Create new folder 
#===========================================
print('CREATE NEW FOLDERS')

# create meta folder, if this doesn't already exist
XenDirectory = "/home/pi/Desktop/XenDirectory/"

if not os.path.exists(XenDirectory):
    os.makedirs(XenDirectory)
os.chdir(XenDirectory)



# New meta folder naming parameters - set above in Configuration 
newMeta = cameraID + audioID
# Create new meta folder if doesn't already exist
print('creating new meta folder')
if os.path.exists(newMeta):
    print('pi + audio folder already exists')
else:
    os.makedirs(newMeta)
    print('pi + audio folder created')



# create observation folder
os.chdir(os.path.join(bettaDirectory,newMeta))
# grab current date
folderDate = dt.datetime.now().strftime('%Y%m%d')
# new observation folder naming parameters
newFolder = "ID" + female + "_pi_" + cameraID + "_audio_"+ audioID + "_" + obs + "_" + folderDate



# create new observation folder
print('creating new folder for {}'.format(obs))
if not os.path.exists(newFolder):
    os.makedirs(newFolder)
    print('new {} folder created'.format(obs))
   

 
# specify file path
filePath = os.path.join(XenDirectory,newMeta,newFolder)



# Update system clock
def update_system_clock():
	cmd = 'sudo systemctl restart systemd-networkd'
	subprocess.call(cmd.split(' '))


update_system_clock()
# Check clock was just updated (should be time very similar to your phone's time)
cmd_check = 'stat /var/lib/systemd/clock'
time_last_update = subprocess.check_output(cmd_check.split(' ')).decode('utf-8').split('\n')[-3:-2]
print('\n')
print('Time updated')
print(time_last_update[0])
print('\n')



#===========================================
# Record video
#===========================================

# Check how long before scheduled start of recording
# -------------------------
todays_date = dt.date.today()

# Schedule start
scheduled_start = dt.datetime.combine(todays_date, dt.time(start_hour,start_min,0,0))
seconds_to_record = (scheduled_start - dt.datetime.now()).seconds
#--------------------------


# If more than a minute before start sleep until 1 minute before start
if seconds_to_record > 60:
	time.sleep(seconds_to_record - 60)

# Within a minute of scheduled start
# Update system clock again
update_system_clock()

# While loop until scheduled start
while dt.datetime.now() < scheduled_start:
	pass


print('RECORDING VIDEO')


#starting video count - leave at zero!
vidCount = 0

omxp = Popen(['omxplayer','-o','local','--loop', movie_path], stdin=subprocess.PIPE,stdout=None, stderr=None, bufsize=0)

with picamera.PiCamera() as camera:
    camera.resolution = (1280,720)
    camera.framerate = 20
    camera.rotation = 270
    camera.start_preview(alpha=5)

    # timestamp settings
    camera.annotate_background = picamera.Color('black')
    camera.annotate_text = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
    
    vidSeg = int(vidHours*60/segLength)

    for i in range (vidSeg):
        vidTime = strftime('%Y%m%d_%H%M%S')
        # create file names for images
        videoIdentifier = "ID" + female + "_pi_" + cameraID + "_audio_"+ audioID + "_" + obs + "_"
        videoName = videoIdentifier + vidTime + ".h264"
        completeVideoFilePath = filePath + videoName

        camera.start_recording(completeVideoFilePath)
        # add start time overlay
        start = dt.datetime.now()
        print(start.strftime('%Y-%m-%d %H:%M:%S'))
        
        while (dt.datetime.now()-start).seconds < segLength*60:
            camera.annotate_text = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
        
        camera.stop_recording()        
       
        # advance vid count - just to keep track of how many videos have been taken
        vidCount += 1
        print('{} videos recorded'.format(vidCount))

    camera.stop_preview()

print('recording complete')

print('GOODBYE')

omxp.stdin.write(b'q')
omxp = Popen(['omxplayer','-i', movie_path])
        



