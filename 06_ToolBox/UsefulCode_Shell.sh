# Shell toolbox

# Change image file names using EXIF date

exiftool -d '%Y%m%d-%H%M-%S%%-02.c.%%e' '-filename<CreateDate' *.JPG
#exiftool -T -FileName -DateTimeOriginal *.JPG > date_time.txt
#exiftool -d '%Y%m%d-%H%M-%S.%%e' '-filename<CreateDate' *.JPG

