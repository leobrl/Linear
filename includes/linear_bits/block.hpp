#pragma once

#include<vector>
#include<memory>

namespace linear{
    
    template <typename T> class Block;
    template <typename T> std::ostream& operator<<( std::ostream&, const Block<T>& );

    template<typename T>
    class Block{
        
        using raw_memory    = std::vector<T>;
        using raw_memory_it = typename std::vector<T>::iterator;

        protected:

            raw_memory buffer;

        public:

            inline                  Block()  = default;
            inline                  ~Block() = default;

            inline explicit         Block(const size_t);
            inline explicit         Block(const std::vector<T>& );

                                    Block(const Block&) = default;
                                    Block(Block&&)      = default;

            Block&                  operator= (const Block&) = default;
            Block&                  operator= (Block&&)      = default;    

            void                    fill(const T&);
            void                    fill(const std::vector<T>&);
            
            inline const T&         operator[] (size_t idx) const   { return buffer[idx]; }
            inline T&               operator[] (size_t idx)         { return buffer[idx]; }

            inline size_t           size () const { return buffer.size(); }

            inline raw_memory_it    begin()    {return buffer.begin();}
            inline raw_memory_it    end()      {return buffer.end();}

            // Friends
            friend std::ostream& operator<< <T> (std::ostream&, const Block<T>&);
    };
}
