import os
import subprocess

compFlag = 0 # 1 for PC, 0 for Mac

if compFlag:
    # PC Path
    vidpath = 'G:\Shared drives\Buhl-Weinburd\Videos\\'
    Hillston = vidpath + 'Hillston Nov 2010\\'
    savepath = 'Data\\examples'
else:
    # Mac Path 
    vidpath = '/Volumes/GoogleDrive/Shared drives/Buhl-Weinburd/Videos/'
    Hillston = vidpath + 'Hillston Nov 2010/'
    savepath = 'Data/examples/'

## 133 video
vid133 = 'VS00133noteLightVariation.MTS'
vid1 = '133_1min_225'
vid2 = '133_1min_650'
vid3 = '133_1min_750'
vid4 = '133_1min_850'
vid5 = '133_1min_950'
vid6 = '133_1min_1115'
vid_Lands = '133_10sec_710_720' # The video that Landsberg tracked manually
vid_Sharma = '133_10sec_1910_1920' # The video that Sharma tracked manually

crop = "crop=w=1824:h=1026:x=48:y=27"
#t0list = ['02:25', '06:50', '7:50', '8:50', '9:50', '11:15']
t0list = ['07:10', '19:10']
#dur = '60' #durations, preceded by -t
dur = '10'

# ### 098 video
# vid098 = 'VS00098.MTS'
# vid1 = '098_1min_430'
# vid2 = '098_1min_530'
# vid3 = '098_1min_630'
# vid4 = '098_1min_730'
# vid5 = '098_1min_5830'
# vid6 = '098_1min_5930'
# vid7 = '098_1min_10030'

# crop = 'crop=1894:976:4:36'
# t0list = ['04:30', '05:30', '06:30', '07:30', '58:30', '59:30', '01:00:30']
# dur = '60'

# ### 096 video
# vid096 = '00096.MTS'
# vid1 = '096_1min_500'
# vid2 = '096_1min_600'
# vid3 = '096_1min_700'
# vid4 = '096_1min_800'
# vid5 = '096_1min_900'
# vid6 = '096_1min_1000'
# vid7 = '096_1min_1100'

# crop = 'crop=1496:802:350:160'
# t0list = ['05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', ]
# dur = '60'

### 146 video
# vid146 = '00146_Low_density.MTS'
# vid1 = '146_1min_515'
# vid2 = '146_1min_615'
# vid3 = '146_1min_715'
# vid4 = '146_1min_815'
# vid5 = '146_1min_915'
# vid6 = '146_1min_1015'
# vid7 = '146_1min_1115'

# crop = 'crop=1650:844:212:30'
# t0list = ['05:15', '06:15', '07:15', '08:15', '09:15', '10:15', '11:15', ]
# dur = '60'

# A list for looping below
# vidlist = [vid1, vid2, vid3, vid4, vid5, vid6, vid7]
vidlist = ['133_10sec_710', '133_10sec_1910']

# For t0list:
# times in the format HH:MM:SS.MILLISECONDS
# can ignore HH, e.g. 00:01 = the first second
# can ignore MM, e.g. 01 = the first second
# can ignore MILLISECONDS
    # or alternatively
    # t0 = #start time
    # tend = #end time, preceded by -to

# For dur:
# duration measured in seconds

# For crop:
# 'crop = width:height:x:y
    # width is of the new output rectangle
    # heights is of the new output rectangle
    # x:y is the new top left corner
    # defaults to (in_w-out_w)/2:(in_h-out_h)/2
# https://ffmpeg.org/ffmpeg-filters.html#crop



# Our video processing, calls ffmpeg
def vidprocess(vidin, t0, dur, crop, vidout):
    '''A routine that calls ffmpeg over and over to preprocess our videos including:
    1)  a) clip from t0 to tend
        b) deinterlace with yadif
        c) encode with libx264 as a .mp4 file
    #2) color inversion
    #3) brightness, contrast, saturation adjustment
    #4) blur
    5) crop to 1824x1026 (these dimensions for the 133 video, based on Landsberg's width=0.95*width)
    #6) greyscale
    7) convert to .AVI format ready for ImageJ to read
    
    Input is a file name without extension, file should be a .MP4 but that may not matter.
    
    Returns Null but creates a series of video files with final version in .AVI format'''

    subprocess.call(['ffmpeg', '-i', vidin, '-ss', t0, '-t', dur, '-vf', 'yadif=0', '-an', '-vcodec',
                    'libx264', '-crf', '10', '-y', 'outfile0.mp4'])
    # -crf 10 gives very good quality. "a subjectively sane range is 17–28"
    #wait = input("Press Enter to continue.") #for debugging
    # invert colors
    #subprocess.call(['ffmpeg','-i', 'outfile0.mp4', '-vf', 'negate', '-y', 'outfile1.mp4'])
    # adjust brightness, contrast
    #subprocess.call(['ffmpeg','-i', 'outfile1.mp4', '-vf', 'eq=brightness=.1:contrast=1000:saturation=1', '-c:a', 'copy', '-y', 'outfile2.mp4' ])
    # blur
    #subprocess.call(['ffmpeg','-i', 'outfile2.mp4', '-vf', "boxblur=5:1", '-y', 'outfile3.mp4'])
    
    # crop
    subprocess.call(['ffmpeg','-i', 'outfile0.mp4', '-vf', crop, '-y', 'outfile4.mp4'])
    # greyscale
    #subprocess.call(['ffmpeg','-i', 'outfile4.mp4', '-vf', 'hue=s=0', '-y', 'outfile5.mp4'])
    
    # reformat
    subprocess.call(['ffmpeg','-i', 'outfile4.mp4', '-pix_fmt', 'nv12', '-f', 'avi', '-vcodec', 'rawvideo', vidout])

    # Delete all the intermediate files
    subprocess.call(['rm','outfile*'])

for [vidin,t0] in zip(vidlist,t0list):
    
    vidfile = Hillston + vid133 # the input video file
    #vidfile = vidin+'.MTS'
    #vidfile = vidin+'.mp4'
    vidout = savepath + 'preprocessed_'+vidin+'.avi' # the output video file
    
    #wait = input("Press Enter to continue.") #for debugging
    vidprocess(vidfile,t0,dur,crop,vidout)
