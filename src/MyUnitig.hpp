#ifndef MYUNITIG_HPP_
#define MYUNITIG_HPP_
#include <ColoredCDBG.hpp>
typedef unsigned int uint32;
class MyUnitig : public CCDBG_Data_t<MyUnitig>, CDBG_Data_t<MyUnitig>
{
public:
    MyUnitig();
    void clear(const UnitigMap<MyUnitig> &um_dest);
    void concat(const UnitigMap<MyUnitig> &um_dest, const UnitigMap<MyUnitig> &um_src);
    void extract(const UnitigMap<MyUnitig> &um_src, bool last_extraction);
    void merge(const UnitigMap<MyUnitig> &um_dest, const const_UnitigMap<MyUnitig> &um_src);
    void clear(const UnitigColorMap<MyUnitig> &um_dest);
    void concat(const UnitigColorMap<MyUnitig> &um_dest, const UnitigColorMap<MyUnitig> &um_src);
    void extract(const UnitigColors *uc_dest, const UnitigColorMap<MyUnitig> &um_src, const bool last_extraction);
    void merge(const UnitigColors &uc_dest, const UnitigColorMap<MyUnitig> &um_dest, const const_UnitigColorMap<MyUnitig> &um_src);
    void toString() const;
    inline void set_visited() { b = 0x1; }
    inline void set_seen() { b = 0x2; }
    inline void set_not_seen_visited() { b = 0x0; }
    inline bool is_visited() const { return (b == 0x1); }
    inline bool is_not_visited() const { return (b != 0x1); }
    inline bool is_seen() const { return (b == 0x2); }
    inline bool is_not_seen() const { return (b != 0x2); }
    inline void set_id(size_t id) { this->id = id; }
    inline size_t get_id() { return id; }
    inline void set_plus_self()
    {
        plus = this;
        set_plus_visited();
    }
    inline void set_minus_self()
    {
        minus = this;
        set_minus_visited();
    }
    inline void set_plus(MyUnitig *p)
    {
        plus = p;
        b = b | 0x01;
    }
    inline void set_minus(MyUnitig *p)
    {
        minus = p;
        b = b | 0x02;
    }
    MyUnitig *get_plus() { return plus; }
    MyUnitig *get_minus() { return minus; }
    MyUnitig *get_plus_or_minus(bool b) { return (b == true) ? plus : minus; }
    inline bool is_plus_super() { return plus != NULL && plus != this; }
    inline bool is_minus_super() { return minus != NULL && minus != this; }
    inline void set_non_super()
    {
        b = (b | 0x04);
    }
    inline bool is_super()
    {
        return !is_both_visited();
    }
    inline bool is_non_super() { return ((0x04) == (b & 0x04)); }
    inline void setComplex(bool strand)
    {
        if (strand == true)
        {
            b = (b | 0x40);
        }
        else
        {
            b = (b | 0x20);
        }
    }
    inline bool is_both_not_null()
    {
        return plus != NULL && minus != NULL;
    }
    inline bool isComplex(bool strand)
    {
        if (strand == true)
        {
            return ((b & 0x40) == (0x40));
        }
        else
        {
            return ((b & 0x20) == (0x20));
        }
    }
    inline bool is_both_visited() { return (b & 0x03) == (0x00); };
    inline bool is_plus_visited() { return ((b & 0x01) == (0x00)); }
    inline bool is_minus_visited() { return ((b & 0x02) == (0x00)); }
    inline void set_plus_visited() { b = (b & (0xFE)); }
    inline void set_minus_visited() { b = (b & (0xFD)); }
    inline void set_both_visited() { b = (b & 0xFC); }
    inline size_t get_bubble_id(bool b)
    {
        return (b == true) ? (*plus).get_id() : (*minus).get_id();
    }
    inline void setStrict(bool strand)
    {
        if (strand == true)
        {
            b = (b | (0x10));
        }
        else
        {
            b = (b | (0x08));
        }
    }
    inline void setNotStrict(bool strand)
    {
        if (strand == true)
        {
            b = (b & (0xEF));
        }
        else
        {
            b = (b & (0xF7));
        }
    }
    inline bool isStrict(bool strand)
    {
        if (strand == true)
        {
            return ((b & 0x10) == (0x10));
        }
        else
        {
            return ((b & 0x08) == (0x08));
        }
    }

private:
    uint8_t b = 0x0;
    MyUnitig *plus = NULL;
    MyUnitig *minus = NULL;
    size_t id;
};
#endif
