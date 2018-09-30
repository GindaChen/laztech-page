# Two_fluid code acceleration

[Blue waters](https://bluewaters.ncsa.illinois.edu/): A [PetaFLOPS-level](https://en.wikipedia.org/wiki/FLOPS) supercomputer located in UIUC. which is fast.

Until 9/30, you will have the permission to log in a series of bluewaters node. It would be great to get our stuff done by the date.

These are the essential questions you will have to answer to get a good start. People are usually struggle by some basic stuff, when they realize that they're actually in the documentation, so it's quite hard to measure how the tutorial should include things and vice versa. 

Instead, try to answer these questions by reading documentations and stuff to set you going. 

## Questions

(Before this step, try to find @Mikie for the necessary access of bluewaters).

1. **Access Machine (fast)**: some basic Shell, `ssh`, `sshpass`, `rsync`

   1. (**Connect the node**) How to login a bluewaters node using `ssh`? How to not enter ssh password every time but instead using [sshpass](https://gist.github.com/arunoda/7790979) to bypass the password? How to use `sshpass` to read a file as the password and `ssh` into the machine? Write this down as a function using Shell, and store them into your `~/.bashrc` file if you're using `bash`.

   2. (**File Transfer**) Usually people use `scp` to upload/download things. But the logic in Bluewaters is different. [This document](https://bluewaters.ncsa.illinois.edu/data-transfer-doc)  doesn't help us because we're just using training accounts, not normal accounts. Read this and see what's the logic behind:

      ```shell
      REMOTEPATH="machine.edu:/my/remote/path"
      LOCALPATH="/my/local/path"
      upload(){
      	rsync -re 'ssh user@machine.edu ssh' $* $REMOTEPATH
      }
      # USAGE: upload a b c d
      
      download(){
      	rsync -re 'ssh user@machine.edu ssh' $REMOTEPATH $LOCALPATH
      }
      # USAGE: download
      ```

      Can you also use `sshpass` to bypass the password for these function? Write these function in a formal way in your `~/.bashrc`.

   3. 

2. **Understand Basic architecture of Supercomputer **

3. *[Important]* **GPU Framework**: OpenACC (OpenMP, maybe)

   1. Read 

4. **Job Submission & Batch Script on Bluewaters**

5. **Programming Language**: Fortran (It's simple)