/**
 * @file
 * Description du fichier
 */

#include <stdbool.h>

/**
 * \brief A salad structure
 */
typedef struct salad
{
  int salad_leaf_count; /** Number of salad leaves */
  int tomato_count; /** Number of tomatoes */
  bool with_vinagre; /** Would you like some vinagre with you salad sir ? */
}
salad_t;

/**
 * Description de la fonction init
 */
int init();

/**
 * Description de la fonction data
 *
 * @param data
 *   Description du parameter data
 *
 * @return
 *   Retourne une valeur X si tout se passe bien
 *   Retourne une valeur différente en cas d'erreur
 *
 * @see processing2()
 */
int processing(data);
