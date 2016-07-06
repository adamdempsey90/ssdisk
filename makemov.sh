ffmpeg -i dens%04d.png -c:v libx264 -vf fps=60 -pix_fmt yuv420p out.mp4
