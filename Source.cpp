#include <iostream>
#include <cstring>
#include <limits> 
#include <climits>
#include <iomanip> // For formatting
#include <thread>  // For sleep simulation
#include <chrono>  // For time
#include <string>  // Explicitly include string

using namespace std;

// ==========================================
//           UI & COLOR UTILITIES
// ==========================================
namespace Color {
    const string RED = "\033[1;31m";
    const string GREEN = "\033[1;32m";
    const string YELLOW = "\033[1;33m";
    const string CYAN = "\033[1;36m";
    const string MAGENTA = "\033[1;35m";
    const string BLUE = "\033[1;34m";
    const string RESET = "\033[0m";
    const string BOLD = "\033[1m";
}

class ConsoleUI {
public:
    static void clearScreen() {
        // ANSI escape code to clear screen
        cout << "\033[2J\033[1;1H";
    }

    static void printHeader(string title) {
        clearScreen();
        cout << Color::CYAN << "========================================================" << Color::RESET << endl;
        int len = (int)title.length();
        int padding = (56 - len) / 2;
        cout << string(padding, ' ') << Color::BOLD << Color::MAGENTA << title << Color::RESET << endl;
        cout << Color::CYAN << "========================================================" << Color::RESET << endl << endl;
    }

    static void printSection(string title) {
        cout << endl << Color::YELLOW << ">>> " << title << " <<<" << Color::RESET << endl;
        cout << Color::YELLOW << "---------------------------------" << Color::RESET << endl;
    }

    static void printSuccess(string msg) {
        cout << Color::GREEN << "[SUCCESS] " << msg << Color::RESET << endl;
    }

    static void printError(string msg) {
        cout << Color::RED << "[ERROR] " << msg << Color::RESET << endl;
    }

    static void loadingBar(string taskName) {
        cout << Color::CYAN << taskName << " " << Color::RESET;
        cout << "[";
        for (int i = 0; i < 20; i++) {
            cout << "="; // Changed to ASCII '=' to avoid encoding warnings
            cout.flush();
            this_thread::sleep_for(chrono::milliseconds(20));
        }
        cout << "] " << Color::GREEN << "DONE!" << Color::RESET << endl;
    }

    static void waitForEnter() {
        cout << endl << Color::CYAN << "[Press Enter to Continue...]" << Color::RESET;
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cin.get();
    }
};

// ==========================================
//           CORE DATA STRUCTURES
// ==========================================

template<typename T1, typename T2>
class Pair {
public:
    T1 first;
    T2 second;

    Pair() : first(T1()), second(T2()) {}
    Pair(const T1& f, const T2& s) : first(f), second(s) {}

    bool operator==(const Pair& other) const {
        return first == other.first && second == other.second;
    }

    bool operator<(const Pair& other) const {
        if (first < other.first) return true;
        if (first > other.first) return false;
        return second < other.second;
    }
};


template<typename T>
class DynamicArray {
private:
    T* data;
    int capacity;
    int length;

    void resize(int newCapacity) {
        T* newData = new T[newCapacity];
        for (int i = 0; i < length; i++) {
            newData[i] = data[i];
        }
        delete[] data;
        data = newData;
        capacity = newCapacity;
    }

public:
    DynamicArray() : data(nullptr), capacity(0), length(0) {}

    DynamicArray(int initialCapacity) : capacity(initialCapacity), length(0) {
        data = new T[capacity];
    }

    ~DynamicArray() {
        if (data) delete[] data;
    }

    DynamicArray(const DynamicArray& other) : capacity(other.capacity), length(other.length) {
        data = new T[capacity];
        for (int i = 0; i < length; i++) {
            data[i] = other.data[i];
        }
    }

    DynamicArray& operator=(const DynamicArray& other) {
        if (this != &other) {
            if (data) delete[] data;
            capacity = other.capacity;
            length = other.length;
            data = new T[capacity];
            for (int i = 0; i < length; i++) {
                data[i] = other.data[i];
            }
        }
        return *this;
    }

    void push_back(const T& value) {
        if (length == capacity) {
            resize(capacity == 0 ? 1 : capacity * 2);
        }
        data[length++] = value;
    }

    void remove_at(int index) {
        for (int i = index; i < length - 1; i++) {
            data[i] = data[i + 1];
        }
        length--;
    }

    T& operator[](int index) { return data[index]; }
    const T& operator[](int index) const { return data[index]; }

    int size() const { return length; }
    bool empty() const { return length == 0; }
    void clear() { length = 0; }
};


template<typename T>
class BaseContainer {
public:
    virtual int size() const = 0;
    virtual bool empty() const = 0;
    virtual void clear_all() = 0;
    virtual void display() const = 0;
    virtual ~BaseContainer() {}
};


template<typename T>
class DiscreteSet : public BaseContainer<T> {
private:
    DynamicArray<T> elements;

    int binary_search(const T& element) const {
        int left = 0;
        int right = elements.size() - 1;
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (elements[mid] == element) return mid;
            else if (elements[mid] < element) left = mid + 1;
            else right = mid - 1;
        }
        return -1;
    }

    int find_insert_position(const T& element) const {
        int left = 0;
        int right = elements.size() - 1;
        int result = elements.size();
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (elements[mid] >= element) {
                result = mid;
                right = mid - 1;
            }
            else {
                left = mid + 1;
            }
        }
        return result;
    }

public:
    DiscreteSet() {}

    void add(const T& element) {
        if (member(element)) return;
        int pos = find_insert_position(element);
        elements.push_back(element);
        for (int i = elements.size() - 1; i > pos; i--) {
            T temp = elements[i];
            elements[i] = elements[i - 1];
            elements[i - 1] = temp;
        }
    }

    void delete_element(const T& element) {
        int index = binary_search(element);
        if (index != -1) elements.remove_at(index);
    }

    bool member(const T& element) const {
        return binary_search(element) != -1;
    }

    int size() const override { return elements.size(); }
    bool empty() const override { return elements.empty(); }
    void clear_all() override { elements.clear(); }
    const T& get_element(int index) const { return elements[index]; }

    DiscreteSet<T> set_union(const DiscreteSet<T>& other) const {
        DiscreteSet<T> result;
        int i = 0, j = 0;
        while (i < this->size() && j < other.size()) {
            if (this->elements[i] < other.elements[j]) result.add(this->elements[i++]);
            else if (this->elements[i] > other.elements[j]) result.add(other.elements[j++]);
            else { result.add(this->elements[i]); i++; j++; }
        }
        while (i < this->size()) result.add(this->elements[i++]);
        while (j < other.size()) result.add(other.elements[j++]);
        return result;
    }

    DiscreteSet<T> set_intersection(const DiscreteSet<T>& other) const {
        DiscreteSet<T> result;
        for (int i = 0; i < this->size(); i++) {
            if (other.member(this->elements[i])) result.add(this->elements[i]);
        }
        return result;
    }

    DiscreteSet<T> set_difference(const DiscreteSet<T>& other) const {
        DiscreteSet<T> result;
        for (int i = 0; i < this->size(); i++) {
            if (!other.member(this->elements[i])) result.add(this->elements[i]);
        }
        return result;
    }

    bool is_subset(const DiscreteSet<T>& other) const {
        for (int i = 0; i < this->size(); i++) {
            if (!other.member(this->elements[i])) return false;
        }
        return true;
    }

    void display() const override {
        cout << Color::CYAN << "{";
        for (int i = 0; i < elements.size(); i++) {
            cout << elements[i];
            if (i < elements.size() - 1) cout << ", ";
        }
        cout << "}" << Color::RESET;
    }
};


class GraphConstraintBuilder {
public:
    static DiscreteSet<int> compute_mandatory_nodes(
        const DiscreteSet<int>& critical_infrastructure,
        const DiscreteSet<int>& high_traffic,
        const DiscreteSet<int>& redundancy_nodes) {

        DiscreteSet<int> union_set = critical_infrastructure.set_union(high_traffic);
        DiscreteSet<int> mandatory = union_set.set_difference(redundancy_nodes);
        return mandatory;
    }
};


template<typename T>
class BinaryRelation {
private:
    DiscreteSet<T> domain;
    DynamicArray<Pair<T, T>> pairs;

public:
    BinaryRelation(const DiscreteSet<T>& dom) : domain(dom) {}

    void add_pair(const T& a, const T& b) {
        Pair<T, T> p(a, b);
        if (!contains_pair(a, b)) pairs.push_back(p);
    }

    bool contains_pair(const T& a, const T& b) const {
        Pair<T, T> p(a, b);
        for (int i = 0; i < pairs.size(); i++) {
            if (pairs[i] == p) return true;
        }
        return false;
    }

    bool is_reflexive() const {
        for (int i = 0; i < domain.size(); i++) {
            if (!contains_pair(domain.get_element(i), domain.get_element(i))) return false;
        }
        return true;
    }

    bool is_symmetric() const {
        for (int i = 0; i < pairs.size(); i++) {
            if (!contains_pair(pairs[i].second, pairs[i].first)) return false;
        }
        return true;
    }

    bool is_transitive() const {
        for (int i = 0; i < pairs.size(); i++) {
            for (int j = 0; j < pairs.size(); j++) {
                if (pairs[i].second == pairs[j].first) {
                    if (!contains_pair(pairs[i].first, pairs[j].second)) return false;
                }
            }
        }
        return true;
    }

    bool is_equivalence() const { return is_reflexive() && is_symmetric() && is_transitive(); }
    bool is_partial_order() const { return is_reflexive() && is_transitive(); }


    bool is_antisymmetric() const {
        for (int i = 0; i < pairs.size(); i++) {
            T a = pairs[i].first;
            T b = pairs[i].second;
            if (a != b && contains_pair(b, a)) return false;
        }
        return true;
    }

    void display() const {
        cout << Color::MAGENTA << "R = {";
        for (int i = 0; i < pairs.size(); i++) {
            cout << "(" << pairs[i].first << ", " << pairs[i].second << ")";
            if (i < pairs.size() - 1) cout << ", ";
        }
        cout << "}" << Color::RESET << endl;
    }
};


template<typename T>
class DiscreteFunction {
private:
    DiscreteSet<T> domain;
    DiscreteSet<T> codomain;
    DynamicArray<Pair<T, T>> mapping;

public:
    DiscreteFunction(const DiscreteSet<T>& dom, const DiscreteSet<T>& codom)
        : domain(dom), codomain(codom) {
    }

    void add_mapping(const T& from, const T& to) {
        Pair<T, T> p(from, to);
        mapping.push_back(p);
    }

    bool is_function() const {
        for (int i = 0; i < domain.size(); i++) {
            T elem = domain.get_element(i);
            int count = 0;
            for (int j = 0; j < mapping.size(); j++) {
                if (mapping[j].first == elem) count++;
            }
            if (count != 1) return false;
        }
        return true;
    }

    bool is_injective() const {
        if (!is_function()) return false;
        for (int i = 0; i < mapping.size(); i++) {
            for (int j = i + 1; j < mapping.size(); j++) {
                if (mapping[i].second == mapping[j].second) return false;
            }
        }
        return true;
    }

    bool is_surjective() const {
        if (!is_function()) return false;
        for (int i = 0; i < codomain.size(); i++) {
            T elem = codomain.get_element(i);
            bool found = false;
            for (int j = 0; j < mapping.size(); j++) {
                if (mapping[j].second == elem) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }

    bool is_bijective() const { return is_injective() && is_surjective(); }

    void display() const {
        cout << "Mapping: " << Color::CYAN << "{";
        for (int i = 0; i < mapping.size(); i++) {
            cout << mapping[i].first << " -> " << mapping[i].second;
            if (i < mapping.size() - 1) cout << ", ";
        }
        cout << "}" << Color::RESET << endl;
    }
};


const int MAX_VERTICES = 100;
void swapInt(int& a, int& b) {
    int temp = a;
    a = b;
    b = temp;
}

struct EdgeNode {
    int destinationVertexId;
    int weight;
    EdgeNode* next;
    EdgeNode(int d = 0, int w = 0) {
        destinationVertexId = d;
        weight = w;
        next = nullptr;
    }
};

struct QNode {
    int key;
    QNode* next;
    QNode(int k = 0) {
        key = k;
        next = nullptr;
    }
};

class Queue {
private:
    QNode* front;
    QNode* rear;
public:
    Queue() { front = nullptr; rear = nullptr; }
    bool isEmpty() { return (front == nullptr); }
    void enqueue(int vid) {
        QNode* n = new QNode(vid);
        if (isEmpty()) { front = rear = n; }
        else { rear->next = n; rear = n; }
    }
    int dequeue() {
        if (isEmpty()) return -1;
        QNode* temp = front;
        int vid = temp->key;
        if (front == rear) front = rear = nullptr;
        else front = front->next;
        delete temp;
        return vid;
    }
    ~Queue() {
        while (front) {
            QNode* n = front;
            front = front->next;
            delete n;
        }
        rear = nullptr;
    }
};

struct HeapNode {
    int vIndex;
    int key;
    HeapNode(int v = -1, int k = INT_MAX) {
        vIndex = v;
        key = k;
    }
};

class MinHeap {
private:
    HeapNode* harr;
    int capacity;
    int heapSize;
    int* pos;
    int posCapacity = 0;

    void swapNode(int i, int j) {
        HeapNode tmp = harr[i];
        harr[i] = harr[j];
        harr[j] = tmp;
        if (harr[i].vIndex >= 0 && harr[i].vIndex < posCapacity) pos[harr[i].vIndex] = i;
        if (harr[j].vIndex >= 0 && harr[j].vIndex < posCapacity) pos[harr[j].vIndex] = j;
    }

    void minHeapify(int i) {
        int l = left(i);
        int r = right(i);
        int smallest = i;
        if (l < heapSize && harr[l].key < harr[smallest].key) smallest = l;
        if (r < heapSize && harr[r].key < harr[smallest].key) smallest = r;
        if (smallest != i) {
            swapNode(i, smallest);
            minHeapify(smallest);
        }
    }

public:
    MinHeap(int cap = MAX_VERTICES * 2, int maxVertices = MAX_VERTICES) {
        capacity = cap;
        harr = new HeapNode[capacity];
        heapSize = 0;
        posCapacity = maxVertices;
        pos = new int[posCapacity];
        for (int i = 0; i < posCapacity; ++i) pos[i] = -1;
    }

    ~MinHeap() {
        delete[] harr;
        delete[] pos;
    }

    bool empty() const { return heapSize == 0; }
    int left(int i) const { return (2 * i + 1); }
    int right(int i) const { return (2 * i + 2); }
    int parent(int i) const { return ((i - 1) / 2); }

    void insertKey(int vIndex, int key) {
        if (heapSize >= capacity) return;
        int i = heapSize++;
        harr[i].vIndex = vIndex;
        harr[i].key = key;
        if (vIndex >= 0 && vIndex < posCapacity) pos[vIndex] = i;
        while (i != 0 && harr[parent(i)].key > harr[i].key) {
            swapNode(i, parent(i));
            i = parent(i);
        }
    }

    void decreaseKey(int vIndex, int newKey) {
        if (vIndex < 0 || vIndex >= posCapacity) return;
        int i = pos[vIndex];
        if (i == -1) return;
        harr[i].key = newKey;
        while (i != 0 && harr[parent(i)].key > harr[i].key) {
            swapNode(i, parent(i));
            i = parent(i);
        }
    }

    HeapNode extractMin() {
        if (heapSize <= 0) return HeapNode(-1, INT_MAX);
        HeapNode root = harr[0];
        if (heapSize == 1) heapSize--;
        else {
            harr[0] = harr[heapSize - 1];
            if (harr[0].vIndex >= 0 && harr[0].vIndex < posCapacity) pos[harr[0].vIndex] = 0;
            heapSize--;
            minHeapify(0);
        }
        if (root.vIndex >= 0 && root.vIndex < posCapacity) pos[root.vIndex] = -1;
        return root;
    }

    bool isInHeap(int vIndex) {
        if (vIndex < 0 || vIndex >= posCapacity) return false;
        return (pos[vIndex] != -1);
    }
};

struct DiscreteSetInt {
    int elements[MAX_VERTICES];
    int size;

    DiscreteSetInt() : size(0) { size = 0; }

    void add(int element) {
        for (int i = 0; i < size; ++i) {
            if (elements[i] == element) return;
        }
        if (size < MAX_VERTICES) {
            elements[size++] = element;
            int i = size - 1;
            while (i > 0 && elements[i - 1] > elements[i]) {
                swapInt(elements[i], elements[i - 1]);
                i--;
            }
        }
    }

    bool member(int element) const {
        int low = 0;
        int high = size - 1;
        while (low <= high) {
            int mid = low + (high - low) / 2;
            if (elements[mid] == element) return true;
            if (elements[mid] < element) low = mid + 1;
            else high = mid - 1;
        }
        return false;
    }

    friend ostream& operator<<(ostream& os, const DiscreteSetInt& s) {
        os << Color::CYAN << "{";
        for (int i = 0; i < s.size; ++i) {
            os << s.elements[i] << (i < s.size - 1 ? "," : "");
        }
        os << "}" << Color::RESET;
        return os;
    }
};

class DSU {
private:
    int parent[MAX_VERTICES];
    int rank[MAX_VERTICES];
    int numVertices;

public:
    DSU(int n = 0) {
        numVertices = n;
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    int find(int i) {
        if (parent[i] == i) return i;
        return (parent[i] = find(parent[i]));
    }

    bool unite(int i, int j) {
        int rooti = find(i);
        int rootj = find(j);
        if (rooti != rootj) {
            if (rank[rooti] < rank[rootj]) parent[rooti] = rootj;
            else if (rank[rooti] > rank[rootj]) parent[rootj] = rooti;
            else {
                parent[rootj] = rooti;
                rank[rooti]++;
            }
            return true;
        }
        return false;
    }

    bool allMandatoryConnected(const DiscreteSetInt& mandatory_vertices) {
        if (mandatory_vertices.size <= 1) return true;
        int first_mandatory_id = -1;
        for (int i = 0; i < mandatory_vertices.size; ++i) {
            int v = mandatory_vertices.elements[i];
            if (v >= 0 && v < numVertices) {
                first_mandatory_id = find(v);
                break;
            }
        }
        if (first_mandatory_id == -1) return true;
        for (int i = 0; i < mandatory_vertices.size; ++i) {
            int v = mandatory_vertices.elements[i];
            if (v >= 0 && v < numVertices) {
                if (find(v) != first_mandatory_id) return false;
            }
        }
        return true;
    }
};

struct KruskalEdge {
    int u, v, weight;
};

void swapKruskalEdge(KruskalEdge* a, KruskalEdge* b) {
    KruskalEdge t = *a;
    *a = *b;
    *b = t;
}

int partition(KruskalEdge arr[], int low, int high) {
    int pivot = arr[high].weight;
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (arr[j].weight <= pivot) {
            i++;
            swapKruskalEdge(&arr[i], &arr[j]);
        }
    }
    swapKruskalEdge(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quicksort(KruskalEdge arr[], int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quicksort(arr, low, pi - 1);
        quicksort(arr, pi + 1, high);
    }
}

struct EdgePair {
    int u, v;
};

struct SpanningTreeResult {
    EdgePair tree_edges[MAX_VERTICES - 1];
    int edgesSize;
    int totalCost;
    bool success;
};

struct ShortestPathResult {
    int vertex_sequence[MAX_VERTICES];
    int sequenceSize;
    int pathCost;
    bool found;
};

class NetworkGraph {
private:
    int vertexCount;
    int edgeCount;
    bool isDirected;
    int matrix[MAX_VERTICES][MAX_VERTICES];
    EdgeNode* lists[MAX_VERTICES];

    void bfs_util(int v, bool visited[]) const;
    bool is_safe_ham(int v, int pos, int path[]) const;
    bool ham_cycle_util(int pos, int path[]) const;

public:
    NetworkGraph(int n, bool directed = false);
    ~NetworkGraph();

    void insert_edge(int from, int to, int weight);
    int edge_weight(int from, int to) const { return matrix[from][to]; }

    int connected_component_count() const;
    bool is_complete_graph() const;
    bool has_euler_circuit() const;
    bool has_euler_path() const;
    bool contains_hamiltonian_cycle() const;
    void generate_property_report() const;

    SpanningTreeResult kruskals_mst(const DiscreteSetInt& mandatory_vertices, int root_vertex) const;
    ShortestPathResult dijkstra_path(int source_vertex, int target_vertex) const;
};

NetworkGraph::NetworkGraph(int n, bool directed) {
    vertexCount = n;
    edgeCount = 0;
    isDirected = directed;
    if (n > MAX_VERTICES) {
        cerr << "Error: Vertex count exceeds MAX_VERTICES. Truncating." << endl;
        vertexCount = MAX_VERTICES;
    }
    for (int i = 0; i < vertexCount; ++i) {
        for (int j = 0; j < vertexCount; ++j) {
            matrix[i][j] = 0;
        }
        lists[i] = nullptr;
    }
}

NetworkGraph::~NetworkGraph() {
    for (int i = 0; i < vertexCount; ++i) {
        EdgeNode* cur = lists[i];
        while (cur) {
            EdgeNode* next = cur->next;
            delete cur;
            cur = next;
        }
    }
}

void NetworkGraph::insert_edge(int from, int to, int weight) {
    if (from < 0 || from >= vertexCount || to < 0 || to >= vertexCount || weight <= 0) return;
    if (matrix[from][to] != 0) {
        matrix[from][to] = weight;
        for (EdgeNode* cur = lists[from]; cur; cur = cur->next) {
            if (cur->destinationVertexId == to) {
                cur->weight = weight;
                break;
            }
        }
        if (!isDirected) {
            matrix[to][from] = weight;
            for (EdgeNode* cur = lists[to]; cur; cur = cur->next) {
                if (cur->destinationVertexId == from) {
                    cur->weight = weight;
                    break;
                }
            }
        }
        return;
    }
    matrix[from][to] = weight;
    EdgeNode* n1 = new EdgeNode(to, weight);
    n1->next = lists[from];
    lists[from] = n1;
    edgeCount++;
    if (!isDirected) {
        matrix[to][from] = weight;
        EdgeNode* n2 = new EdgeNode(from, weight);
        n2->next = lists[to];
        lists[to] = n2;
    }
}

void NetworkGraph::bfs_util(int v, bool visited[]) const {
    Queue q;
    visited[v] = true;
    q.enqueue(v);
    while (!q.isEmpty()) {
        int u = q.dequeue();
        for (EdgeNode* cur = lists[u]; cur; cur = cur->next) {
            int neighbor = cur->destinationVertexId;
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.enqueue(neighbor);
            }
        }
    }
}

int NetworkGraph::connected_component_count() const {
    if (vertexCount == 0) return 0;
    bool visited[MAX_VERTICES];
    for (int i = 0; i < vertexCount; ++i) visited[i] = false;
    int count = 0;
    for (int i = 0; i < vertexCount; ++i) {
        bool has_edge = false;
        for (int j = 0; j < vertexCount; ++j) {
            if (matrix[i][j] != 0 || matrix[j][i] != 0) {
                has_edge = true;
                break;
            }
        }
        if (has_edge && !visited[i]) {
            bfs_util(i, visited);
            count++;
        }
        else if (!has_edge && !visited[i]) {
            count++;
            visited[i] = true;
        }
    }
    return count;
}

bool NetworkGraph::is_complete_graph() const {
    if (vertexCount <= 1) return true;
    int expected_edges = isDirected ? vertexCount * (vertexCount - 1) : vertexCount * (vertexCount - 1) / 2;
    int current_edges = isDirected ? edgeCount : edgeCount / 2;
    return current_edges == expected_edges;
}

bool NetworkGraph::has_euler_circuit() const {
    if (vertexCount == 0 || connected_component_count() > 1) return false;
    if (!isDirected) {
        for (int i = 0; i < vertexCount; ++i) {
            int degree = 0;
            for (EdgeNode* cur = lists[i]; cur; cur = cur->next) degree++;
            if (degree % 2 != 0) return false;
        }
        return true;
    }
    return false;
}

bool NetworkGraph::has_euler_path() const {
    if (has_euler_circuit()) return true;
    if (vertexCount == 0 || connected_component_count() > 1) return false;
    if (!isDirected) {
        int odd_degree_count = 0;
        for (int i = 0; i < vertexCount; ++i) {
            int degree = 0;
            for (EdgeNode* cur = lists[i]; cur; cur = cur->next) degree++;
            if (degree % 2 != 0) odd_degree_count++;
        }
        return (odd_degree_count == 2);
    }
    return false;
}

bool NetworkGraph::is_safe_ham(int v, int pos, int path[]) const {
    if (matrix[path[pos - 1]][v] == 0) return false;
    for (int i = 0; i < pos; ++i) {
        if (path[i] == v) return false;
    }
    return true;
}

bool NetworkGraph::ham_cycle_util(int pos, int path[]) const {
    if (pos == vertexCount) return (matrix[path[pos - 1]][path[0]] != 0);
    for (int v = 0; v < vertexCount; ++v) {
        if (is_safe_ham(v, pos, path)) {
            path[pos] = v;
            if (ham_cycle_util(pos + 1, path)) return true;
            path[pos] = -1;
        }
    }
    return false;
}

bool NetworkGraph::contains_hamiltonian_cycle() const {
    if (vertexCount < 3) return false;
    int path[MAX_VERTICES];
    for (int i = 0; i < vertexCount; ++i) path[i] = -1;
    path[0] = 0;
    return ham_cycle_util(1, path);
}

void NetworkGraph::generate_property_report() const {
    cout << Color::MAGENTA << "\n[ NETWORK GRAPH PROPERTIES ]" << Color::RESET << endl;
    cout << "---------------------------------" << endl;
    cout << left << setw(25) << "Total Vertices:" << vertexCount << endl;
    cout << left << setw(25) << "Total Edges:" << (isDirected ? edgeCount : edgeCount / 2) << endl;
    cout << left << setw(25) << "Graph Type:" << (isDirected ? "[DIRECTED]" : "[UNDIRECTED]") << endl;
    cout << left << setw(25) << "Euler Circuit:" << (has_euler_circuit() ? Color::GREEN + "TRUE" + Color::RESET : Color::RED + "FALSE" + Color::RESET) << endl;
    cout << left << setw(25) << "Euler Path:" << (has_euler_path() ? Color::GREEN + "TRUE" + Color::RESET : Color::RED + "FALSE" + Color::RESET) << endl;
    cout << left << setw(25) << "Hamiltonian Cycle:" << (contains_hamiltonian_cycle() ? Color::GREEN + "TRUE" + Color::RESET : Color::RED + "FALSE" + Color::RESET) << endl;
    cout << left << setw(25) << "Complete Graph:" << (is_complete_graph() ? Color::GREEN + "TRUE" + Color::RESET : Color::RED + "FALSE" + Color::RESET) << endl;
    cout << left << setw(25) << "Connected Components:" << connected_component_count() << endl;
}

SpanningTreeResult NetworkGraph::kruskals_mst(const DiscreteSetInt& mandatory_vertices, int root_vertex) const {

    SpanningTreeResult result;
    result.edgesSize = 0;
    result.totalCost = 0;
    result.success = false;

    if (vertexCount == 0) {
        result.success = true;
        return result;
    }

    KruskalEdge edges[MAX_VERTICES * MAX_VERTICES];
    int edge_list_size = 0;

    for (int i = 0; i < vertexCount; ++i) {
        int start_j = isDirected ? 0 : i + 1;
        for (int j = start_j; j < vertexCount; ++j) {
            if (matrix[i][j] != 0) {
                edges[edge_list_size++] = { i, j, matrix[i][j] };
            }
        }
    }

    quicksort(edges, 0, edge_list_size - 1);
    DSU dsu(vertexCount);
    int edges_in_mst = 0;

    for (int i = 0; i < edge_list_size; ++i) {
        int u = edges[i].u;
        int v = edges[i].v;
        int w = edges[i].weight;
        if (dsu.unite(u, v)) {
            result.tree_edges[result.edgesSize++] = { u, v };
            result.totalCost += w;
            edges_in_mst++;
        }
    }

    if (vertexCount > 0 && edges_in_mst != vertexCount - 1) {
        result.success = false;
        return result;
    }
    if (!dsu.allMandatoryConnected(mandatory_vertices)) {
        result.success = false;
        return result;
    }
    result.success = true;
    return result;
}

ShortestPathResult NetworkGraph::dijkstra_path(int source_vertex, int target_vertex) const {

    ShortestPathResult result;
    result.sequenceSize = 0;
    result.pathCost = INT_MAX;
    result.found = false;

    if (source_vertex < 0 || source_vertex >= vertexCount || target_vertex < 0 || target_vertex >= vertexCount) return result;
    if (source_vertex == target_vertex) {
        result.vertex_sequence[0] = source_vertex;
        result.sequenceSize = 1;
        result.pathCost = 0;
        result.found = true;
        return result;
    }

    int dist[MAX_VERTICES];
    int parent[MAX_VERTICES];
    MinHeap mh(MAX_VERTICES * 2, MAX_VERTICES);

    for (int i = 0; i < vertexCount; ++i) {
        dist[i] = INT_MAX;
        parent[i] = -1;
        if (i < vertexCount) mh.insertKey(i, INT_MAX);
    }

    dist[source_vertex] = 0;
    mh.decreaseKey(source_vertex, 0);

    while (!mh.empty()) {
        HeapNode hn = mh.extractMin();
        int u = hn.vIndex;
        int du = hn.key;
        if (du == INT_MAX) break;
        if (u == target_vertex) break;

        for (EdgeNode* cur = lists[u]; cur; cur = cur->next) {
            int v = cur->destinationVertexId;
            int weight = cur->weight;
            if (mh.isInHeap(v)) {
                if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    parent[v] = u;
                    mh.decreaseKey(v, dist[v]);
                }
            }
        }
    }

    if (dist[target_vertex] != INT_MAX) {
        result.found = true;
        result.pathCost = dist[target_vertex];
        int temp_path[MAX_VERTICES];
        int curr = target_vertex;
        int temp_size = 0;
        while (curr != -1) {
            temp_path[temp_size++] = curr;
            curr = parent[curr];
        }
        for (int i = 0; i < temp_size; ++i) {
            result.vertex_sequence[i] = temp_path[temp_size - 1 - i];
        }
        result.sequenceSize = temp_size;
    }
    return result;
}



class IntegrationController {
private:
    void clear_cin() {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    void parse_set_from_line(const char* line, DiscreteSet<int>& set) {
        int num = 0;
        bool in_number = false;
        bool is_negative = false;
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == '-' && !in_number) {
                is_negative = true;
                in_number = true;
            }
            else if (line[i] >= '0' && line[i] <= '9') {
                num = num * 10 + (line[i] - '0');
                in_number = true;
            }
            else if (in_number) {
                if (is_negative) num = -num;
                set.add(num);
                num = 0;
                in_number = false;
                is_negative = false;
            }
        }
        if (in_number) {
            if (is_negative) num = -num;
            set.add(num);
        }
    }

    DiscreteSet<int> get_user_set(const char* name) {
        cout << "Enter elements for " << Color::YELLOW << name << Color::RESET << " (space separated, e.g., 1 2 3): ";
        char buffer[1024];
        cin.getline(buffer, 1024);
        DiscreteSet<int> s;
        parse_set_from_line(buffer, s);
        return s;
    }

public:
    void run_graph_constraint_module() {
        ConsoleUI::printHeader("GRAPH CONSTRAINT BUILDER");

        DiscreteSet<int> critical = get_user_set("Critical Infrastructure Nodes");
        DiscreteSet<int> high_traffic = get_user_set("High Traffic Nodes");
        DiscreteSet<int> redundancy = get_user_set("Redundancy Nodes");

        cout << endl;
        ConsoleUI::loadingBar("Processing Constraints");

        DiscreteSet<int> mandatory = GraphConstraintBuilder::compute_mandatory_nodes(
            critical, high_traffic, redundancy);

        cout << Color::GREEN << "Mandatory Nodes: " << Color::RESET;
        mandatory.display();
        cout << endl;
        ConsoleUI::waitForEnter();
    }

    void run_set_operations_module() {
        ConsoleUI::printHeader("SET OPERATIONS ANALYZER");

        DiscreteSet<int> A = get_user_set("Set A");
        DiscreteSet<int> B = get_user_set("Set B");

        ConsoleUI::printSection("Results");
        cout << left << setw(15) << "Union:" << Color::RESET; A.set_union(B).display(); cout << endl;
        cout << left << setw(15) << "Intersection:" << Color::RESET; A.set_intersection(B).display(); cout << endl;
        cout << left << setw(15) << "Difference:" << Color::RESET; A.set_difference(B).display(); cout << endl;

        ConsoleUI::waitForEnter();
    }

    void run_relations_module() {
        ConsoleUI::printHeader("BINARY RELATION ANALYZER");

        DiscreteSet<int> domain = get_user_set("Relation Domain");
        BinaryRelation<int> rel(domain);

        cout << endl << Color::CYAN << "Enter pairs (enter 'x' or non-integer to stop entry):" << Color::RESET << endl;
        while (true) {
            int a, b;
            cout << "Pair " << Color::YELLOW << "(a b)" << Color::RESET << ": ";
            if (!(cin >> a >> b)) { clear_cin(); break; }
            rel.add_pair(a, b);
        }

        cout << endl;
        ConsoleUI::loadingBar("Analyzing Properties");
        rel.display();

        // FIXED: Added '+' for string concatenation
        cout << left << setw(15) << "Reflexive:" << (rel.is_reflexive() ? Color::GREEN + "TRUE" : Color::RED + "FALSE") << Color::RESET << endl;
        cout << left << setw(15) << "Symmetric:" << (rel.is_symmetric() ? Color::GREEN + "TRUE" : Color::RED + "FALSE") << Color::RESET << endl;
        cout << left << setw(15) << "Transitive:" << (rel.is_transitive() ? Color::GREEN + "TRUE" : Color::RED + "FALSE") << Color::RESET << endl;

        ConsoleUI::waitForEnter();
    }

    void run_function_module() {
        ConsoleUI::printHeader("FUNCTION PROPERTIES ANALYZER");

        DiscreteSet<int> domain = get_user_set("Domain");
        DiscreteSet<int> codomain = get_user_set("Codomain");
        DiscreteFunction<int> func(domain, codomain);

        cout << endl << Color::CYAN << "Enter mappings (enter 'x' or non-integer to stop entry):" << Color::RESET << endl;
        while (true) {
            int a, b;
            cout << "Map " << Color::YELLOW << "(from to)" << Color::RESET << ": ";
            if (!(cin >> a >> b)) { clear_cin(); break; }
            func.add_mapping(a, b);
        }

        cout << endl;
        ConsoleUI::loadingBar("Checking Function Rules");
        func.display();

        bool isFunc = func.is_function();
        // FIXED: Added '+' for string concatenation
        cout << left << setw(15) << "Is Function:" << (isFunc ? Color::GREEN + "TRUE" : Color::RED + "FALSE") << Color::RESET << endl;
        if (isFunc) {
            cout << left << setw(15) << "Is Bijective:" << (func.is_bijective() ? Color::GREEN + "TRUE" : Color::RED + "FALSE") << Color::RESET << endl;
        }

        ConsoleUI::waitForEnter();
    }


    void run_interactive_graph_analysis() {
        ConsoleUI::printHeader("CDAF GRAPH ANALYZER");

        int numVertices;
        cout << "Enter the total number of vertices (e.g., 5): ";
        if (!(cin >> numVertices)) { clear_cin(); return; }

        NetworkGraph g(numVertices, false);

        cout << endl << Color::CYAN << ">>> INPUT EDGES <<<" << Color::RESET << endl;
        cout << "Format: " << Color::YELLOW << "Source Destination Weight" << Color::RESET << endl;
        cout << "Type " << Color::RED << "-1 -1 -1" << Color::RESET << " to finish adding edges." << endl;

        int u, v, w;
        while (true) {
            cout << "Edge: ";
            cin >> u >> v >> w;
            if (u == -1) break;
            g.insert_edge(u, v, w);
        }

        ConsoleUI::loadingBar("Generating Property Report");
        g.generate_property_report();

        ConsoleUI::printSection("Kruskal's MST Calculation");
        DiscreteSetInt mandatory_nodes;
        cout << "Enter mandatory nodes for MST (enter -1 to stop): ";
        int node;
        while (cin >> node && node != -1) {
            mandatory_nodes.add(node);
        }

        SpanningTreeResult mst_result = g.kruskals_mst(mandatory_nodes, 0);
        cout << "Mandatory Vertices: " << mandatory_nodes << endl;

        // Sort for display consistency
        for (int i = 0; i < mst_result.edgesSize; ++i) {
            for (int j = i + 1; j < mst_result.edgesSize; ++j) {
                int u1 = mst_result.tree_edges[i].u;
                int v1 = mst_result.tree_edges[i].v;
                int u2 = mst_result.tree_edges[j].u;
                int v2 = mst_result.tree_edges[j].v;
                if (u1 > v1) swapInt(u1, v1);
                if (u2 > v2) swapInt(u2, v2);
                if (u2 < u1 || (u2 == u1 && v2 < v1)) {
                    EdgePair temp = mst_result.tree_edges[i];
                    mst_result.tree_edges[i] = mst_result.tree_edges[j];
                    mst_result.tree_edges[j] = temp;
                }
            }
        }

        if (mst_result.success) {
            for (int i = 0; i < mst_result.edgesSize; ++i) {
                int u = mst_result.tree_edges[i].u;
                int v = mst_result.tree_edges[i].v;
                if (u > v) swapInt(u, v);
                cout << "Edge (" << u << ", " << v << "): Weight " << Color::GREEN << g.edge_weight(u, v) << Color::RESET << endl;
            }
            cout << Color::BOLD << "TOTAL COST: " << mst_result.totalCost << Color::RESET << endl;
        }
        else {
            ConsoleUI::printError("MST Construction Failed (Connectivity or Constraints issue)");
        }

        ConsoleUI::printSection("Dijkstra's Shortest Path");
        int source, target;
        cout << "Enter Source and Target vertices (e.g., 0 4): ";
        cin >> source >> target;

        ShortestPathResult path_result = g.dijkstra_path(source, target);

        cout << "SHORTEST PATH (" << source << " -> " << target << "): ";
        if (path_result.found) {
            cout << Color::GREEN << "[ ";
            for (int i = 0; i < path_result.sequenceSize; ++i) {
                cout << path_result.vertex_sequence[i] << (i < path_result.sequenceSize - 1 ? " -> " : "");
            }
            cout << " ]" << Color::RESET << endl;
            cout << "TOTAL COST: " << Color::BOLD << path_result.pathCost << Color::RESET << endl;
        }
        else {
            ConsoleUI::printError("NO PATH FOUND");
            cout << "TOTAL COST: INFINITY" << endl;
        }
        clear_cin();
        ConsoleUI::waitForEnter();
    }

    void start_menu() {
        while (true) {
            ConsoleUI::clearScreen();
            cout << Color::CYAN << R"(
   ______   ___       ____    ______     _______.____    ____  _______.
  /      | |   \     /    \  |   ____|   /       |\   \  /   / /       |
 |  ,----' |    \   /  ^   \ |  |__     |   (----` \   \/   / |   (----`
 |  |      |  .  \ /  /_\   \|   __|     \   \      \_    _/   \   \    
 |  `----. |  |\  \  /  _   \  |  |  .----)   |       |  | .----)   |   
  \______| |__| \__\/__/ \__\__|  |__|_______/        |__| |_______/    
            )" << Color::RESET << endl;

            cout << Color::YELLOW << "   Discrete Structures & Graph Analysis Integration System" << Color::RESET << endl;
            cout << Color::BLUE << "   Developed by [Usman Ali]  :) " << Color::RESET << endl;
            cout << endl;
            cout << "================================================================" << endl;
            cout << "  " << Color::BOLD << "[1]" << Color::RESET << " Set Operations (Union, Intersection, Diff)" << endl;
            cout << "  " << Color::BOLD << "[2]" << Color::RESET << " Graph Constraints (Mandatory Nodes)" << endl;
            cout << "  " << Color::BOLD << "[3]" << Color::RESET << " Binary Relations (Properties Analysis)" << endl;
            cout << "  " << Color::BOLD << "[4]" << Color::RESET << " Function Properties (Injection, Surjection)" << endl;
            cout << "  " << Color::BOLD << "[5]" << Color::RESET << " Advanced Graph Analysis (MST, Dijkstra, Euler)" << endl;
            cout << "  " << Color::BOLD << "[0]" << Color::RESET << " Exit System" << endl;
            cout << "================================================================" << endl;
            cout << Color::GREEN << "  Select Option: " << Color::RESET;

            int choice;
            if (!(cin >> choice)) {
                clear_cin();
                continue;
            }
            clear_cin();

            if (choice == 0) {
                cout << Color::CYAN << "\nExiting... Goodbye!" << Color::RESET << endl;
                break;
            }

            switch (choice) {
            case 1: run_set_operations_module(); break;
            case 2: run_graph_constraint_module(); break;
            case 3: run_relations_module(); break;
            case 4: run_function_module(); break;
            case 5: run_interactive_graph_analysis(); break;
            default: ConsoleUI::printError("Invalid Option Selected."); this_thread::sleep_for(chrono::seconds(1));
            }
        }
    }
};

int main() {
    IntegrationController controller;
    controller.start_menu();
    return 0;
}