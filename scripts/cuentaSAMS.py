import csv
import os
from collections import Counter
from itertools import combinations
import math

def calcular_jaccard(set1, set2):
    """Calcula el índice de Jaccard entre dos conjuntos"""
    interseccion = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return interseccion / union if union > 0 else 0

def calcular_shannon(frecuencias):
    """Calcula el índice de diversidad de Shannon"""
    total = sum(frecuencias.values())
    shannon = 0
    for count in frecuencias.values():
        if count > 0:
            p = count / total
            shannon -= p * math.log(p)
    return shannon

def analizar_heterogeneidad(archivo_csv):
    """Analiza la heterogeneidad transcriptómica del archivo CSV"""
    
    if not os.path.exists(archivo_csv):
        print(f"Error: El archivo '{archivo_csv}' no existe.\n")
        return None
    
    try:
        # Estructuras de datos
        sams_por_paciente = {}  # {paciente_id: set(SAMs)}
        pacientes_por_sam = {}  # {SAM: set(pacientes)}
        distribucion_frecuencias = Counter()  # {n_pacientes: count}
        
        # Leer el archivo CSV
        with open(archivo_csv, 'r', encoding='utf-8') as file:
            csv_reader = csv.DictReader(file)
            
            for row in csv_reader:
                locus = row['locus']
                n_patients = int(row['n_patients'])
                patients_str = row['patients']
                
                # Parsear la lista de pacientes
                pacientes = [p.strip() for p in patients_str.split(',')]
                
                # Actualizar estructuras
                pacientes_por_sam[locus] = set(pacientes)
                distribucion_frecuencias[n_patients] += 1
                
                for paciente in pacientes:
                    if paciente not in sams_por_paciente:
                        sams_por_paciente[paciente] = set()
                    sams_por_paciente[paciente].add(locus)
        
        # 1. ÍNDICE DE JACCARD (promedio entre todos los pares de pacientes)
        print("\n" + "="*70)
        print("CÁLCULO DEL ÍNDICE DE JACCARD")
        print("="*70)
        
        pacientes_lista = list(sams_por_paciente.keys())
        num_pacientes = len(pacientes_lista)
        print(f"Número total de pacientes: {num_pacientes}")
        
        if num_pacientes < 2:
            print("Se necesitan al menos 2 pacientes para calcular Jaccard")
            jaccard_promedio = None
        else:
            jaccard_scores = []
            # Calcular Jaccard para todos los pares de pacientes
            for p1, p2 in combinations(pacientes_lista, 2):
                jaccard = calcular_jaccard(sams_por_paciente[p1], sams_por_paciente[p2])
                jaccard_scores.append(jaccard)
            
            jaccard_promedio = sum(jaccard_scores) / len(jaccard_scores)
            jaccard_min = min(jaccard_scores)
            jaccard_max = max(jaccard_scores)
            
            print(f"Número de pares comparados: {len(jaccard_scores):,}")
            print(f"\nÍndice de Jaccard promedio: {jaccard_promedio:.4f}")
            print(f"Jaccard mínimo: {jaccard_min:.4f}")
            print(f"Jaccard máximo: {jaccard_max:.4f}")
            print(f"\nInterpretación:")
            print(f"  - Valores cercanos a 0: Alta heterogeneidad (pacientes muy diferentes)")
            print(f"  - Valores cercanos a 1: Baja heterogeneidad (pacientes muy similares)")
        
        # 2. ÍNDICE DE SHANNON
        print("\n" + "="*70)
        print("ÍNDICE DE DIVERSIDAD DE SHANNON")
        print("="*70)
        
        shannon = calcular_shannon(distribucion_frecuencias)
        print(f"Índice de Shannon: {shannon:.4f}")
        print(f"\nInterpretación:")
        print(f"  - Valores altos: Mayor diversidad en la distribución de SAMs")
        print(f"  - Valores bajos: SAMs concentrados en pocas frecuencias")
        
        # 3. DISTRIBUCIÓN DE FRECUENCIAS
        print("\n" + "="*70)
        print("DISTRIBUCIÓN DE FRECUENCIAS")
        print("="*70)
        print(f"{'Nº Pacientes':<15} {'Nº SAMs':<15} {'Porcentaje':<15}")
        print("-"*70)
        
        total_sams = sum(distribucion_frecuencias.values())
        for n_pac in sorted(distribucion_frecuencias.keys()):
            count = distribucion_frecuencias[n_pac]
            porcentaje = (count / total_sams) * 100
            print(f"{n_pac:<15} {count:<15,} {porcentaje:>6.2f}%")
        
        print("-"*70)
        print(f"{'TOTAL':<15} {total_sams:<15,} {100.0:>6.2f}%")
        
        # Estadísticas adicionales
        print("\n" + "="*70)
        print("ESTADÍSTICAS ADICIONALES")
        print("="*70)
        print(f"Total de SAMs únicos: {len(pacientes_por_sam):,}")
        print(f"Total de pacientes: {num_pacientes}")
        print(f"SAMs por paciente (promedio): {sum(len(sams) for sams in sams_por_paciente.values()) / num_pacientes:.2f}")
        
        # SAMs únicos (solo en 1 paciente)
        sams_unicos = sum(1 for sams in pacientes_por_sam.values() if len(sams) == 1)
        print(f"SAMs únicos (solo 1 paciente): {sams_unicos:,} ({sams_unicos/total_sams*100:.2f}%)")
        
        # SAMs compartidos (en todos los pacientes)
        sams_compartidos = sum(1 for sams in pacientes_por_sam.values() if len(sams) == num_pacientes)
        print(f"SAMs en todos los pacientes: {sams_compartidos:,} ({sams_compartidos/total_sams*100:.2f}%)")
        
        return {
            'jaccard_promedio': jaccard_promedio,
            'shannon': shannon,
            'distribucion': distribucion_frecuencias,
            'total_sams': total_sams,
            'num_pacientes': num_pacientes
        }
    
    except Exception as e:
        print(f"Error al procesar '{archivo_csv}': {e}\n")
        import traceback
        traceback.print_exc()
        return None

def guardar_resultados(archivo_csv, resultados):
    """Guarda los resultados en un archivo de texto"""
    nombre_salida = archivo_csv.replace('.csv', '_heterogeneidad.txt')
    
    with open(nombre_salida, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write(f"ANÁLISIS DE HETEROGENEIDAD TRANSCRIPTÓMICA\n")
        f.write(f"Archivo: {archivo_csv}\n")
        f.write("="*70 + "\n\n")
        
        # Jaccard
        f.write("ÍNDICE DE JACCARD\n")
        f.write("-"*70 + "\n")
        if resultados['jaccard_promedio'] is not None:
            f.write(f"Jaccard promedio: {resultados['jaccard_promedio']:.4f}\n\n")
        else:
            f.write("No calculado (menos de 2 pacientes)\n\n")
        
        # Shannon
        f.write("ÍNDICE DE SHANNON\n")
        f.write("-"*70 + "\n")
        f.write(f"Shannon: {resultados['shannon']:.4f}\n\n")
        
        # Distribución
        f.write("DISTRIBUCIÓN DE FRECUENCIAS\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Nº Pacientes':<15} {'Nº SAMs':<15} {'Porcentaje':<15}\n")
        f.write("-"*70 + "\n")
        
        total_sams = resultados['total_sams']
        for n_pac in sorted(resultados['distribucion'].keys()):
            count = resultados['distribucion'][n_pac]
            porcentaje = (count / total_sams) * 100
            f.write(f"{n_pac:<15} {count:<15} {porcentaje:>6.2f}%\n")
        
        f.write("-"*70 + "\n")
        f.write(f"{'TOTAL':<15} {total_sams:<15} {100.0:>6.2f}%\n")
    
    print(f"\n✓ Resultados guardados en '{nombre_salida}'")

def main():
    print("\n*** ANÁLISIS DE HETEROGENEIDAD TRANSCRIPTÓMICA ***\n")
    
    while True:
        archivo = input("Introduce el nombre del archivo CSV (o 'salir' para terminar): ").strip()
        
        if archivo.lower() == 'salir':
            print("\n¡Análisis finalizado!")
            break
        
        if not archivo.endswith('.csv'):
            archivo += '.csv'
        
        resultados = analizar_heterogeneidad(archivo)
        
        if resultados:
            guardar = input("\n¿Guardar resultados en archivo? (s/n): ").strip().lower()
            if guardar == 's':
                guardar_resultados(archivo, resultados)
            print("\n")

if __name__ == "__main__":
    main()