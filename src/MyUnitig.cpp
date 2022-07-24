#include "MyUnitig.hpp"
#include <queue>
#include <stack>
MyUnitig::MyUnitig() : b(0x00) {} // Initiate the boolean to "not visited"
void MyUnitig::clear(const UnitigMap<MyUnitig> &um_dest)
{
	set_not_seen_visited(); // Set the unitig to "not visited"
}
void MyUnitig::concat(const UnitigMap<MyUnitig> &um_dest, const UnitigMap<MyUnitig> &um_src)
{
	set_not_seen_visited(); // Set the unitig to "not visited"
}
void MyUnitig::extract(const UnitigMap<MyUnitig> &um_src, bool last_extraction)
{
	set_not_seen_visited(); // Set the unitig to "not visited"
}
void MyUnitig::merge(const UnitigMap<MyUnitig> &um_dest, const const_UnitigMap<MyUnitig> &um_src)
{
	set_not_seen_visited(); // Set the unitig to "not visited"
}
void MyUnitig::clear(const UnitigColorMap<MyUnitig> &um_dest)
{
	set_not_seen_visited(); // Set the unitig to "not visited"
}
void MyUnitig::concat(const UnitigColorMap<MyUnitig> &um_dest, const UnitigColorMap<MyUnitig> &um_src)
{
	set_not_seen_visited(); // Set the unitig to "not visited"
}
void MyUnitig::extract(const UnitigColors *uc_dest, const UnitigColorMap<MyUnitig> &um_src, const bool last_extraction)
{
	set_not_seen_visited();
}
void MyUnitig::merge(const UnitigColors &uc_dest, const UnitigColorMap<MyUnitig> &um_dest, const const_UnitigColorMap<MyUnitig> &um_src)
{
	set_not_seen_visited();
}
void MyUnitig::toString() const
{
	cout << "Unitig visited = " << (is_visited() ? "true" : "false") << endl;
	cout << "Unitig seen = " << (is_seen() ? "true" : "false") << endl;
}
