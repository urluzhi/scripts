// assume Unix, use fcntl
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

// file locking (so multiple users don't interfere with each other)
// bWrite = true for write/exclusive access, else read/shared access
// returns a file descriptor if successful, else -1
// note: locks are released when the files are closed at end of program
int LockPlease(const char * filename,const bool bWrite)
{
    const int MAX_TRY = 10;
    int fd; // file descriptor
    struct flock lck;

    fd = open(filename, (bWrite ? O_WRONLY : O_RDONLY));
    if (fd < 0)
        return -1; // file not found, or could not be opened

    // set up the record locking structure, the address of which
    // is passed to the fcntl system call.
    lck.l_type = (bWrite ? F_WRLCK : F_RDLCK);  // write or read lock
    lck.l_whence = 0;                           // from beginning of file
    lck.l_start = lck.l_len = 0L;               // until end of the file

    // Attempt locking MAX_TRY times before giving up.
    for (int tryno = 0; tryno < MAX_TRY; tryno++)
        if (fcntl(fd, F_SETLK, &lck) == 0)
            return fd; // lock obtained OK
        else if (errno == EAGAIN || errno == EACCES) // file in use
            sleep(2); // wait 2 seconds and try again
        else
            return -1; // fcntl error other than file in use
    return -1; // MAX_TRY exceeded
}



//lock the file for share reading(0) in cgi
int lockdat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, char *hexaloop, 
	  char *tstacki23, char *tstacki1n,bool RW) {
	//lock the files using the C++ method for reading share
 LockPlease(miscloop,RW);
 LockPlease(loop2,RW);
 LockPlease(stackf,RW);
 LockPlease(tstackh,RW);
 LockPlease(tstacki,RW);
 LockPlease(tloop,RW);
 LockPlease(danglef,RW);
 LockPlease(int22,RW);
 LockPlease(int21,RW);
 LockPlease(triloop,RW);
 LockPlease(coax,RW);
 LockPlease(tstackcoax,RW);
 LockPlease(coaxstack,RW);
 LockPlease(tstack,RW);
 LockPlease(tstackm,RW);
 LockPlease(int11,RW);
 LockPlease(hexaloop,RW);
 LockPlease(tstacki23,RW);
 LockPlease(tstacki1n,RW);


}

