import numpy as np
import matplotlib.pyplot as plt

#Reads a square image in 8-bit/color PPM format from the given file. Note: No checks on valid format are done.
def readImage(filename):
    f = open(filename,"rb")
    
    f.readline()
    s = f.readline()
    f.readline()
    (pixel, pixel) = [t(s) for t,s in zip((int,int),s.split())]
    
    data = np.fromfile(f,dtype=np.uint8,count = pixel*pixel*3)
    img = data.reshape((pixel,pixel,3)).astype(np.double)
    
    f.close()
    
    return img, pixel
    

#Writes a square image in 8-bit/color PPM format.
def writeImage(filename, image):
    f = open(filename,"wb")
    
    pixel = image.shape[0]
    ppm_header = f'P6 ' + str(pixel) + ' ' + str(pixel) + ' ' + str(255) +'\n'
    f.write(bytearray(ppm_header,'utf-8'))
    
    image = image.astype(np.uint8)
    
    image.tofile(f)
    
    f.close()
    
    
# kernel function in real space
def W_real(r,h):

    k = 40/(7*np.pi*h**2)
    if 0 <= r/h < 1/2:
        return k*(1 - 6*(r/h)**2 + 6*(r/h)**3)
    elif 1/2 <= r/h < 1:
        return k*(2*(1-r/h)**3)
    else:
        return 0
    

img, pixel = readImage("aq-original.ppm")

color_sums_before = np.zeros(3)
for i in range(3):
    color_sums_before[i] = np.sum(img[:,:,i])
print(f"Red pixels sum: {color_sums_before[0]}") 
print(f"Green pixels sum: {color_sums_before[1]}") 
print(f"Blue pixels sum: {color_sums_before[2]}") 


#Now we set up our desired smoothing kernel. We'll use complex number for it even though it is real. 
kernel_real = np.zeros((pixel,pixel),dtype=np.complex)

hsml = 15.

#now set the values of the kernel 
for i in np.arange(pixel):
    for j in np.arange(pixel):
        r = np.sqrt((i-pixel/2)**2 + (j-pixel/2)**2)
        kernel_real[i,j] = W_real(r,hsml)

kernel_real = np.fft.fftshift(kernel_real)        
        #TODO: do something sensible here to set the real part of the kernel
        #kernel_real[i, j] = ....
        


#Let's calculate the Fourier transform of the kernel
kernel_kspace = np.fft.fft2(kernel_real)


#further space allocations for image transforms
color_real = np.zeros((pixel,pixel),dtype=np.complex)


#we now convolve each color channel with the kernel using FFTs
for colindex in np.arange(3):
    #copy input color into complex array
    color_real[:,:].real = img[:,:,colindex]
    
    
    #forward transform
    color_kspace = np.fft.fft2(color_real)
    
    #multiply with kernel in Fourier space
    #TODO: fill in code here

    color_kspace *= kernel_kspace
    
    #backward transform
    color_real = np.fft.ifft2(color_kspace)
    
    
    #copy real value of complex result back into color array
    img[:,:,colindex] = color_real.real
    

color_sums_after = np.zeros(3)
for i in range(3):
    color_sums_after[i] = np.sum(img[:,:,i])
print(f"Red pixels sum: {color_sums_after[0]}") 
print(f"Green pixels sum: {color_sums_after[1]}") 
print(f"Blue pixels sum: {color_sums_after[2]}") 

writeImage("aq-smoothed.ppm", img)
