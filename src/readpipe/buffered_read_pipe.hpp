// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef buffered_read_pipe_hpp
#define buffered_read_pipe_hpp

namespace octopus {

class BufferedReadPipe
{
public:
    BufferedReadPipe() = delete;
    
    BufferedReadPipe(const BufferedReadPipe&)            = delete;
    BufferedReadPipe& operator=(const BufferedReadPipe&) = delete;
    BufferedReadPipe(BufferedReadPipe&&)                 = default;
    BufferedReadPipe& operator=(BufferedReadPipe&&)      = default;
    
    ~BufferedReadPipe() = default;
    
private:

};

} // namespace octopus

#endif
