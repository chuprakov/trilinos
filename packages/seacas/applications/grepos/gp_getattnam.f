C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 
      subroutine getattnam(ndb, nelblk, id, numatr, names)
      include 'gp_namlen.blk'
      integer ndb, nelblk
      integer id(*)
      integer numatr(*)
      character*(maxnam) names(*)
      
      ibeg = 1
      do iel = 1, nelblk
        if (numatr(iel) .gt. 0) then
          call exgean(ndb, id(iel), numatr(iel), names(ibeg), ierr)
        end if
        ibeg = ibeg + numatr(iel)
      end do

      return 
      end

      subroutine putattnam(ndb, nelblk, id, numatr, names)
      include 'gp_namlen.blk'
      integer ndb, nelblk
      integer id(*)
      integer numatr(*)
      character*(maxnam) names(*)
      
      ibeg = 1
      do iel = 1, nelblk
        if (numatr(iel) .gt. 0) then
          call expean(ndb, id(iel), numatr(iel), names(ibeg), ierr)
        end if
        ibeg = ibeg + numatr(iel)
      end do

      return 
      end
